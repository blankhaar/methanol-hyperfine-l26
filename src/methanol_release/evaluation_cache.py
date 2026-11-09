"""Optional exact-content memoization for one original DirectProblem evaluation.

No model, data, weights, arithmetic, derivatives, or step sizes are rewritten.
Original functions run on every distinct input. Cache values are deep-copied on
return, preserving the original functions' fresh-output semantics. In particular,
Dataset.profile_prediction normalizes transition strengths in place.

This context temporarily replaces process-local function bindings. It must only
be used in a single-threaded worker, never concurrently by multiple Python
threads. Separate optimizer processes are safe. Original bindings are restored
on exit, including on exceptions. No project files are edited. Profile dicts
must not be mutated concurrently during an evaluation (also an original-model
requirement). Memory is bounded by a single evaluate/jacobian call.
"""
from contextlib import contextmanager
from functools import lru_cache
import copy
import numpy as np


def array_key(a):
    a=np.asarray(a)
    return (a.dtype.str,a.shape,a.tobytes(order='C'))


def value_key(v):
    if isinstance(v,np.ndarray):return ('array',array_key(v))
    if isinstance(v,(float,np.floating)):return ('float',array_key(np.asarray(v)))
    if isinstance(v,(list,tuple)):return tuple(value_key(q) for q in v)
    return v


def freeze_signature(args,kwargs):
    return value_key(args),tuple((k,value_key(v)) for k,v in sorted(kwargs.items()))


@contextmanager
def exact_evaluation_cache(problem,stats=None,features=None):
    """Cache exact repeated work for the enclosed original model evaluation.

    features optionally selects any subset of gaussian, diagonalize, wigner,
    transitions, components, direct_prediction, low_prediction. All by default.
    stats receives calls/hits/misses/entries for each enabled operation.
    """
    import methanol_hyperfine.spin as spin
    import methanol_hyperfine.lambdip as lambdip
    import methanol_release.forward as dataset
    import methanol_release.forward as direct
    if stats is None:stats={}
    if features is None:features={'gaussian','diagonalize','wigner','transitions','components','direct_prediction','low_prediction'}
    changes=[]
    wigner_cache=None
    def memo(name,fn,key=freeze_signature):
        cache={};count=stats[name]={'calls':0,'hits':0,'misses':0,'entries':0}
        def wrapped(*args,**kwargs):
            count['calls']+=1;k=key(args,kwargs)
            if k in cache:
                count['hits']+=1
                return copy.deepcopy(cache[k])
            value=fn(*args,**kwargs);cache[k]=copy.deepcopy(value)
            count['misses']+=1;count['entries']=len(cache)
            return value
        return wrapped
    def patch(obj,name,fn):
        own=name in obj.__dict__;old=getattr(obj,name);changes.append((obj,name,old,own));setattr(obj,name,fn)
    try:
        if 'gaussian' in features:patch(lambdip,'gaussian_2f',memo('gaussian',lambdip.gaussian_2f))
        if 'wigner' in features:
            # This six-integer public function has scalar immutable output.
            # Use the C-backed LRU key path rather than generic array freezing.
            fn=lru_cache(maxsize=None,typed=True)(spin.wigner3j);wigner_cache=fn
            patch(spin,'wigner3j',fn);patch(lambdip,'wig3jj',fn)
        if 'diagonalize' in features:
            fn=memo('diagonalize',spin.diagonalize)
            patch(spin,'diagonalize',fn);patch(lambdip,'diagonalize',fn);patch(dataset,'diagonalize',fn)
        if 'transitions' in features:
            fn=memo('transitions',spin.transition_components)
            patch(spin,'transition_components',fn);patch(dataset,'transition_components',fn)
        if 'components' in features:
            def components_key(args,kwargs):
                g,h=args
                return g['u'],g['l'],g['species'],array_key(h[g['u']]),array_key(h[g['l']])
            fn=memo('components',lambdip.components,components_key)
            patch(lambdip,'components',fn);patch(direct,'components',fn)
        if 'direct_prediction' in features:
            def prediction_key(args,kwargs):
                g,*rest=args
                return id(g),freeze_signature(rest,kwargs)
            fn=memo('direct_prediction',lambdip.prediction,prediction_key)
            patch(lambdip,'prediction',fn);patch(direct,'prediction',fn)
        if 'low_prediction' in features:
            def low_key(args,kwargs):
                g,h,nuisance=args
                return id(g),array_key(h[g['u']]),array_key(h[g['l']]),array_key(nuisance)
            patch(problem.d,'profile_prediction',memo('low_prediction',problem.d.profile_prediction,low_key))
        yield stats
    finally:
        if wigner_cache is not None:
            info=wigner_cache.cache_info();stats['wigner']={'calls':info.hits+info.misses,'hits':info.hits,'misses':info.misses,'entries':info.currsize}
        for obj,name,old,own in reversed(changes):
            if own:setattr(obj,name,old)
            else:delattr(obj,name)


def install_exact_cache(problem,features=None):
    """Opt one existing problem into exact memoization, retaining its identity.

    Usage: install_exact_cache(problem). All existing residual, jacobian,
    evaluate, and explicit problem._cache=None operations continue to work.
    Reversible with uninstall_exact_cache(problem). Do not install twice.
    """
    if '_exact_cache_installation' in problem.__dict__:
        raise RuntimeError('Exact cache is already installed on this problem')
    original=problem.evaluate;owned='evaluate' in problem.__dict__
    def evaluate(z,jac=False):
        stats={};problem._exact_cache_last_stats=stats
        with exact_evaluation_cache(problem,stats,features):
            return original(z,jac)
    problem._exact_cache_installation=(original,owned)
    problem.evaluate=evaluate
    return problem


def uninstall_exact_cache(problem):
    """Restore exactly the original evaluation binding."""
    original,owned=problem.__dict__.pop('_exact_cache_installation')
    if owned:problem.evaluate=original
    else:del problem.evaluate
    problem.__dict__.pop('_exact_cache_last_stats',None)
    return problem
