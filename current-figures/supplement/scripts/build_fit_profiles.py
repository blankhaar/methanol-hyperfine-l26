from pathlib import Path
import json,glob,hashlib,shutil,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1]
archive=BASE/'data/profile_metadata.json'
meta={'profiles':json.load(open(archive))}
aa=dict(np.load(BASE/'data/fit_profiles.npz'))
endpoint=BASE.parent/'inputs/fit_endpoint.json'
assert hashlib.sha256(endpoint.read_bytes()).hexdigest()==json.loads((BASE/'data/endpoint_pin.json').read_text())['sha256']
res=np.load(BASE/'data/fit_residual.npz')['residual']; delta=float(np.max(abs(res-aa['fit_full_weighted_residual'])))
assert delta<1e-10
# The endpoint is included verbatim in the package inputs/.
profiles={k:v for k,v in meta['profiles'].items()}
json.dump(profiles,open(BASE/'data/profile_metadata.json','w'),indent=2)
arr={k:v for k,v in aa.items() if k.startswith(('direct_','low_','fit_','Belov_all57_'))}
np.savez_compressed(BASE/'data/fit_profiles.npz',**arr)
main={'direct_00','direct_05','low_01','low_03','low_05','low_06','low_10','low_15'}
plt.rcParams.update({'font.size':9,'axes.labelsize':8,'axes.titlesize':9,'font.family':'DejaVu Sans','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
manifest=[]
for family,n,batch in [('direct',18,6),('low',16,4)]:
 for page,start in enumerate(range(0,n,batch),1):
  keys=[f'{family}_{i:02}' for i in range(start,min(start+batch,n))]
  rows=len(keys)//2
  fig=plt.figure(figsize=(7.1,rows*2.6)); gs=fig.add_gridspec(rows,2,hspace=.52,wspace=.32)
  for pos,key in enumerate(keys):
   r,c=divmod(pos,2); inner=gs[r,c].subgridspec(2,1,height_ratios=[3,1],hspace=.06); ax=fig.add_subplot(inner[0]); rr=fig.add_subplot(inner[1],sharex=ax)
   m=profiles[key]; x=aa[key+'_x']; y=aa[key+'_y']; pred=aa['fit_'+key+'_prediction']; mask=aa[key+'_mask'].astype(bool) if key+'_mask' in aa else np.ones(len(x),bool)
   ax.plot(x,y,color='.55',lw=.8,label='Digitized observation');ax.plot(x,pred,color='#087e8b',lw=1.15,label='L26');rr.axhline(0,color='.65',lw=.6);rr.plot(x,y-pred,color='#087e8b',lw=.8)
   if family=='low':
    for l,h in m.get('artifact_intervals_khz',[]):
     ax.axvspan(l,h,color='#e9c46a',alpha=.3);rr.axvspan(l,h,color='#e9c46a',alpha=.3)
    title=f"{m['lab']}: {m['reference_mhz']/1000:.4f} GHz {m['species']}"
   else:title=f"Xu Fig. {m['figure']}: {m['panel']}, K={m['ku']}←{m['kl']}"
   ax.set_title(title+(' *' if key in main else ''),loc='left',pad=5)
   ax.text(.99,.96,m['id'],ha='right',va='top',transform=ax.transAxes,fontsize=6.5,color='.35')
   ax.set_ylabel('Signal');rr.set_ylabel('O−C');rr.set_xlabel('Frequency offset (kHz)');ax.tick_params(labelbottom=False)
   ax.margins(x=.015);rr.margins(x=.015)
   if pos==0:ax.legend(fontsize=6.6,loc='lower left',frameon=False)
   out=np.column_stack([x,y,pred,y-pred,mask.astype(int)])
   np.savetxt(BASE/'data'/f'{m["id"]}.csv',out,delimiter=',',header='frequency_offset_kHz,observation,fit,observed_minus_calculated,included_mask',comments='')
   manifest.append({'key':key,'scan_id':m['id'],'family':family,'figure_page':page,'main_paper_display':key in main,'samples':len(x),'included_samples':int(mask.sum()),'rms_normalized_signal':float(np.sqrt(np.mean((y[mask]-pred[mask])**2)))})
  fig.subplots_adjust(left=.075,right=.965,bottom=.07,top=.965)
  f=BASE/'figures'/f'{"xu" if family=="direct" else "low"}-{page}.pdf';fig.savefig(f,bbox_inches="tight",pad_inches=.05);plt.close(fig)
json.dump({'endpoint_sha256':hashlib.sha256(endpoint.read_bytes()).hexdigest(),'residual_replay_max_abs':delta,'profiles':manifest,'recomputed_fit':False},open(BASE/'data/profile_validation.json','w'),indent=2)
with open(BASE/'data/profile_index.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=manifest[0]);w.writeheader();w.writerows(manifest)
print('Generated 7 atlas plates, 34 profile CSVs; residual replay max',delta)
