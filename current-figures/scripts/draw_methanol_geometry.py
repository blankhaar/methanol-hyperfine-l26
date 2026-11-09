#!/usr/bin/env python3
"""Rebuild the vector methanol structure figure (Python + NumPy + Matplotlib).

Geometry: Lankhaar, Groenenboom & van der Avoird, JCP 145, 244301
(2016), Sec. III B.1, DOI 10.1063/1.4972004. The reported rigid-top
bond lengths and angles are used; the torsional pose is illustrative.

The torsional zero is the staggered configuration with H1 anti to H4.
Positive alpha is the negative right-hand methyl rotation about C->O,
consistent with p_alpha=-i d/dalpha and F(p_alpha+rho P_z)^2.
The diagram defines the BELGI kinetic coordinate. Tensor conventions are
specified independently in docs/MODEL.md.
All primitives are vector paths/text; no raster molecular rendering.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/methanol-mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Arc
from datetime import datetime, timezone
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'manuscript' / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix',
    'font.size': 11, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'svg.fonttype': 'none', 'svg.hashsalt': 'methanol-geometry-belgi-v1',
    'axes.linewidth': .6, 'savefig.facecolor':'white',
})

# Cartesian frame used only to construct the illustration: C->O is +x,
# H4 lies in x-y. No identification with RAM or principal a,b,c axes.
rco, rch, roh = 1.427, 1.096, .956
coh, och = np.deg2rad([108.87, 109.91])
# Preserve the molecular pose; assign its signed BELGI coordinate explicitly.
alpha_deg = 35.
torsional_reference = np.pi
pose = torsional_reference - np.deg2rad(alpha_deg)
positions = {'C': np.array([0.,0.,0.]), 'O': np.array([rco,0.,0.]),
             'H4': np.array([rco-roh*np.cos(coh), roh*np.sin(coh),0.])}
for i in range(1,4):
    theta = pose - (i-1)*2*np.pi/3
    positions[f'H{i}'] = rch*np.array([np.cos(och),np.sin(och)*np.cos(theta),np.sin(och)*np.sin(theta)])

# Diagonalize the full molecular inertia tensor about its center of mass.
# These exact eigenvectors, with a right-handed sign convention, are projected
# using the same camera as every atom and bond.
masses={'C':12.,'O':15.99491461957,**{f'H{i}':1.00782503223 for i in range(1,5)}}
com=sum(masses[k]*v for k,v in positions.items())/sum(masses.values())
I=sum(masses[k]*((q@q)*np.eye(3)-np.outer(q,q)) for k,p in positions.items() for q in [p-com])
moments,U=np.linalg.eigh(I)
a=U[:,0]*np.sign(U[0,0]);b=U[:,1]*np.sign(U[1,1]);c=np.cross(a,b)
PAM=np.array([a,b,c]).T
# Right-handed RAM axes: z parallel to I^-1 lambda, x in the C-O-H4
# plane toward H4, and y=z cross x. This is a geometric construction;
# it neither replaces the fitted rho nor rotates the stored tensor arrays.
lambda_CO=np.array([1.,0.,0.])
z_ram=np.linalg.solve(I,lambda_CO);z_ram/=np.linalg.norm(z_ram)
x_ram=b-(b@z_ram)*z_ram;x_ram/=np.linalg.norm(x_ram)
y_ram=np.cross(z_ram,x_ram)
RAM=np.column_stack([x_ram,y_ram,z_ram])
assert np.max(abs(RAM.T@RAM-np.eye(3)))<1e-12
assert np.linalg.det(RAM)>0
assert z_ram@lambda_CO>0 and x_ram@(positions['H4']-positions['O'])>0
# Oblique orthographic camera, specified in the principal frame. Screen up
# is the projected principal a direction, so a is vertical while b and c
# visibly recede to opposite sides. All atoms, bonds and axes share this view.
view_principal=np.array([.4,.55,.73]);view_principal/=np.linalg.norm(view_principal)
view=PAM@view_principal
screen_y=a-(a@view)*view;screen_y/=np.linalg.norm(screen_y)
screen_x=np.cross(screen_y,view)
rot=np.eye(2)
def project(x):
    x=np.asarray(x)-com
    return np.stack([x@screen_x,x@screen_y],axis=-1)

# Compact vertical stack. Cartesian geometry and the camera are unchanged;
# the displayed axes are now explicitly the right-handed RAM axes.
fig=plt.figure(figsize=(2.7,4.6))
ax=fig.add_axes([.025,.495,.95,.495])
ax.set_aspect('equal'); ax.set_xlim(-1.55,1.55);ax.set_ylim(-1.85,2.10);ax.axis('off')
fig.text(.11,.975,'(a)',ha='left',va='top',fontsize=11)

# Fine axis extension: deliberately not a principal-inertia or RAM axis.
line=project(np.array([[-.95,0,0],[2.70,0,0]]))
ax.plot(line[:,0],line[:,1],color='.6',lw=.7,ls=(0,(3,3)),zorder=0)
la,lb=project(np.array([[2.22,0,0],[2.70,0,0]]))
ax.annotate('',xy=lb,xytext=la,arrowprops={'arrowstyle':'-|>','color':'.35','lw':.8},zorder=1)
ax.text(lb[0]+.13,lb[1]-.05,r'$\boldsymbol{\lambda}$',ha='left',va='center',fontsize=13)

# RAM directions originate at the molecular center of mass. The lengths
# only control the visible arrow extents; their 3D directions and orthographic
# projection are exactly those obtained from this same molecular geometry.
origin=project(com)
axis_lengths={'z':1.90,'x':1.50,'y':1.70}
axis_labels={
    'z':{'offset':np.array([-.12,.045]),'ha':'right','va':'bottom'},
    'x':{'offset':np.array([-.10,-.015]),'ha':'right','va':'center'},
    'y':{'offset':np.array([.10,-.015]),'ha':'left','va':'center'},
}
axis_endpoints={}
for name,vec in zip('zxy',[z_ram,x_ram,y_ram]):
    endpoint=project(com+axis_lengths[name]*vec)
    axis_endpoints[name]=endpoint
    ax.annotate('',xy=endpoint,xytext=origin,
                arrowprops={'arrowstyle':'-|>','color':'.25','lw':.85,
                            'linestyle':(0,(4,2.5)), 'shrinkA':0,'shrinkB':0},
                zorder=3)
    settings=axis_labels[name]
    ax.text(*(endpoint+settings['offset']),rf'${name}$',fontsize=13,
            ha=settings['ha'],va=settings['va'],zorder=20,
            bbox={'facecolor':'white','edgecolor':'none','pad':.25})

# Sphere radii are pictorial, not atomic/van der Waals radii.
radii={'C':.255,'O':.235,'H1':.178,'H2':.178,'H3':.178,'H4':.178}
basegray={'C':.55,'O':.20,'H1':.93,'H2':.93,'H3':.93,'H4':.93}

def sphere(axis,center,radius,base,label,fs=12,z=10):
    center=np.asarray(center)
    # Nested offset circles give restrained gray illumination while remaining
    # native vector objects (no image or mesh rasterization).
    axis.add_patch(Circle(center,radius,facecolor=str(base),edgecolor='none',zorder=z))
    for t in np.linspace(0,1,42)[1:]:
        rr=radius*(1-.74*t)
        offset=radius*np.array([-.17,.19])*t
        gray=base+(1-base)*.55*t
        axis.add_patch(Circle(center+offset,rr,facecolor=str(gray),edgecolor='none',zorder=z+.01*t))
    axis.add_patch(Circle(center,radius,facecolor='none',edgecolor='black',lw=.9,zorder=z+.1))
    axis.text(*center,label,ha='center',va='center',fontsize=fs,color='white' if base<.3 else 'black',zorder=z+.2)

def bond(axis,a,b,z):
    pa,pb=project(a),project(b)
    axis.plot([pa[0],pb[0]],[pa[1],pb[1]],color='black',lw=6.9,solid_capstyle='round',zorder=z)
    axis.plot([pa[0],pb[0]],[pa[1],pb[1]],color='.83',lw=5.15,solid_capstyle='round',zorder=z+.01)
    axis.plot([pa[0],pb[0]],[pa[1],pb[1]],color='.97',lw=1.25,solid_capstyle='round',zorder=z+.02)

# All five projected bonds are disjoint except at atom centers. Draw entire
# sticks beneath opaque atom disks; this avoids spurious depth-segment seams.
for u,v in [('C','O'),('C','H1'),('C','H2'),('C','H3'),('O','H4')]:
    bond(ax,positions[u],positions[v],2)
for k,p in sorted(positions.items(),key=lambda kv:kv[1]@view):
    label=k if k in ['C','O'] else rf'$\mathrm{{H}}_{{{k[-1]}}}$'
    sphere(ax,project(p),radii[k],basegray[k],label,fs=13 if k in ['C','O'] else 12,z=10+p@view)

ax.plot(*origin,marker='o',ms=3.0,mfc='white',mec='.25',mew=.7,zorder=20)

# Group captions and unobtrusive leaders.
ax.text(0.,-1.77,r'$\mathrm{CH_3}$ top',ha='center',va='center',fontsize=11)
ax.text(-1.00,1.20,'OH frame',ha='center',va='center',fontsize=11)
ax.plot([-.98,-.78],[1.07,.95],color='.4',lw=.65,zorder=20)

# Axial inset is essential: perspective does not preserve a dihedral angle.
ins=fig.add_axes([.025,.015,.95,.455]);ins.set_aspect('equal')
ins.set_xlim(-1.40,1.40);ins.set_ylim(-1.75,1.69);ins.axis('off')
fig.text(.11,.465,'(b)',ha='left',va='top',fontsize=11)
ins.text(.13,1.65,r'View along C$\,\rightarrow\,$O',ha='center',va='top',fontsize=10)
# Looking along +x: +y projects upwards and +z to the right.
def radial(theta,r=1.):return r*np.array([np.sin(theta),np.cos(theta)])
# A light outer circle identifies the rear O, the inner C is in front.
ins.add_patch(Circle((0,0),.21,fc='white',ec='black',lw=.8,zorder=2))
# OH radial projection; the central section is hidden by the front atom.
h4=radial(0,1.02)
ins.plot([0,h4[0]],[0,h4[1]],color='black',lw=1.3,zorder=1)
for i in range(1,4):
    h=radial(pose-(i-1)*2*np.pi/3,1.0)
    ins.plot([0,h[0]],[0,h[1]],color='black',lw=1.3,zorder=3)
    ins.add_patch(Circle(h,.105,fc='white',ec='black',lw=.8,zorder=4))
    labelpoint=radial(pose-(i-1)*2*np.pi/3,1.22)
    ins.text(*labelpoint,rf'$\mathrm{{H}}_{{{i}}}$',ha='center',va='center',fontsize=11,zorder=5)
ins.add_patch(Circle(h4,.105,fc='white',ec='black',lw=.8,zorder=4))
ins.text(0,1.20,r'$\mathrm{H}_4$',ha='center',va='center',fontsize=11,zorder=5)
ins.add_patch(Circle((0,0),.07,fc='black',ec='black',lw=.6,zorder=6))
# Zero is staggered: the H1 reference ray is opposite the H4 ray. It is
# a coordinate reference, not an additional atom. Positive alpha is
# counterclockwise in this C-front/O-behind view.
refstart,refend=radial(torsional_reference,.25),radial(torsional_reference,1.02)
ins.plot([refstart[0],refend[0]],[refstart[1],refend[1]],
         color='.4',lw=.85,ls=(0,(3,2.5)),zorder=1)
ins.text(0,-1.23,r'$\mathrm{H}_1(0)$',ha='center',va='center',fontsize=10,color='.25')
angles=np.linspace(torsional_reference,pose,100)
arc=np.array([radial(t,.59) for t in angles])
ins.plot(arc[:,0],arc[:,1],color='black',lw=1.05,zorder=7)
ins.annotate('',xy=arc[-1],xytext=arc[-15],
             arrowprops={'arrowstyle':'-|>','mutation_scale':10,'lw':1.05,
                         'color':'black','shrinkA':0,'shrinkB':0},zorder=8)
al=radial((torsional_reference+pose)/2,.78)
ins.text(*al,r'$\alpha$',fontsize=13,ha='center',va='center',zorder=8,
         bbox={'facecolor':'white','edgecolor':'none','pad':.5})

# At the intended 3.3-inch column width, the smallest label is 8.8 pt.
# Scale type separately from the unchanged molecular construction.
for text in fig.findobj(matplotlib.text.Text):
    text.set_fontsize(.8 * text.get_fontsize())

metadata={'Title':'Methanol: RAM axes and signed BELGI torsional coordinate',
          'Author':'', 'Subject':'Right-handed RAM axes; alpha increases counterclockwise viewed C to O',
          'Keywords':'methanol, hyperfine, torsion, vector illustration',
          'CreationDate':datetime(2026,10,14,tzinfo=timezone.utc),
          'ModDate':datetime(2026,11,8,tzinfo=timezone.utc)}
fig.savefig(OUT/'methanol_geometry.pdf',metadata=metadata)
fig.savefig(OUT/'methanol_geometry.svg',metadata={'Date':'2026-11-08','Creator':'draw_methanol_geometry.py'})
fig.savefig(OUT/'methanol_geometry.png',dpi=300)
plt.close(fig)

provenance={
  'figure':'methanol_geometry.pdf', 'source_script':'scripts/draw_methanol_geometry.py',
  'primary_source':{'citation':'B. Lankhaar, G. C. Groenenboom, and A. van der Avoird, J. Chem. Phys. 145, 244301 (2016), Sec. III B.1',
    'doi':'10.1063/1.4972004','url':'https://publications.lib.chalmers.se/records/fulltext/247517/local_247517.pdf',
    'printed_page':'244301-5','pdf_page_one_based':6},
  'bond_lengths_A':{'CO':rco,'CH':rch,'OH':roh},'bond_angles_deg':{'COH':108.87,'OCH':109.91},
  'illustrative_H1_azimuth_from_OH_plane_deg':145, 'torsional_reference_H1_azimuth_deg':180, 'illustrative_signed_BELGI_alpha_deg':alpha_deg,
  'proton_numbering':{'methyl':[1,2,3],'hydroxyl':4},
  'geometry_scope':'Reported C3-symmetrized rigid-top bond parameters; illustrative non-equilibrium torsional pose; atom radii schematic.',
  'convention_scope':'BELGI plus kinetic convention: p_alpha=-i d/dalpha and F(p_alpha+rho P_z)^2. Lambda points C to O; alpha increases counterclockwise in the C-front/O-behind view and equals minus the right-hand geometric gamma. This figure clarifies the kinetic coordinate; fixed tensor source mapping is separately qualified and no fitted constants or numerical inputs are changed.',
  'proton_azimuth_formula':'theta_i=pi-alpha-2*pi*(i-1)/3, i=1,2,3; H4 azimuth=0; native alpha_i=alpha+2*pi*(i-1)/3',
  'proton_label_revision':'H2 and H3 exchanged relative to the prior drawing to match native positive Fourier phases; molecular positions as a set unchanged.',
  'construction_coordinates_A':{k:v.tolist() for k,v in positions.items()},
  'checks':{'methyl_bond_lengths':all(abs(np.linalg.norm(positions[f'H{i}'])-rch)<1e-12 for i in range(1,4)),
            'OH_length':abs(np.linalg.norm(positions['H4']-positions['O'])-roh)<1e-12,
            'CO_length':abs(np.linalg.norm(positions['O'])-rco)<1e-12},
  'raster_content_in_pdf':False,
  'axes_display':'RAM x,y,z originate at the molecular center of mass; z parallel to I^-1 lambda, x in C-O-H4 plane toward H4, y=z cross x. Actual directions use the same camera as the atoms; arrow lengths are schematic. Principal axes are retained only for construction/provenance, not displayed or identified with RAM.',
  'RAM_axes_columns_xyz':RAM.tolist(),
  'RAM_z_to_PAM_a_deg':float(np.rad2deg(np.arccos(np.clip(z_ram@a,-1,1)))),
  'RAM_z_to_CO_deg':float(np.rad2deg(np.arccos(np.clip(z_ram@lambda_CO,-1,1)))),
  'axis_arrow_lengths_A':axis_lengths, 'axis_origin_projected':origin.tolist(),
  'axis_endpoints_projected':{k:v.tolist() for k,v in axis_endpoints.items()},
  'projection':{'view_toward_observer':view.tolist(),'view_in_principal_frame':view_principal.tolist(),'screen_x':screen_x.tolist(),'screen_y':screen_y.tolist(),'world_origin':com.tolist(),'in_plane_rotation_deg':0,'screen_up':'exact projected principal a direction'},
  'masses_u':masses, 'center_of_mass_A':com.tolist(), 'inertia_uA2':I.tolist(),
  'principal_moments_uA2':moments.tolist(), 'principal_axes_columns':PAM.tolist(),
  'lambda_display':'Unit direction from C to O; z is parallel to rho=I_alpha I^-1 lambda.',
  'kinetic_convention_source':{'citation':'Jansen et al., Phys. Rev. A 84, 062505 (2011), Eqs. (5d), (13)', 'url':'https://arxiv.org/html/1109.5076v1'},
  'fixed_tensor_qualification':'The L26 uses the tested DD pullback D_ij(alpha)=D_source_pi(i)pi(j)(-alpha), pi=(2 3), and corrected single-c sine terms of M. Native BELGI wavefunctions and p_alpha are unchanged. The drawing itself is independent of fitted constants.',
  'torsional_zero_primary_evidence':[{'source':'Lankhaar et al. JCP145,244301 (2016)','locations':'SecIII B1 equilibrium HOCH=180deg; Eq25 positive V3 potential; SecII C (23)* reflection fixes H1 and H4','doi':'10.1063/1.4972004'},{'source':'Coudert et al. Phys Rev A100,043425 (2019)','locations':'SecIII A, p043425-3 and Fig3: explicit gamma=0 staggered H1 antiH4','doi':'10.1103/PhysRevA.100.043425','url':'https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevA.100.043425/fulltext'}],
}
(ROOT/'evidence'/'methanol_geometry_provenance.json').write_text(json.dumps(provenance,indent=2,default=lambda obj:obj.item() if isinstance(obj,np.generic) else str(obj))+'\n')
print('Created',OUT/'methanol_geometry.pdf')
