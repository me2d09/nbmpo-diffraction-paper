import numpy as np
import xml.etree.ElementTree as ET
from fractions import Fraction

from dataclasses import dataclass
import matplotlib.pyplot as plt

# Data files, fill your paths correctly
PATH_XML = '.'
PATH_OVF = rf'???'
PATH_VESTA = rf'???'


results = [
(0, r'\mag_600mK_0T\Na2BaMnP2O8-600mK-0T-Hc-ball-ppp.xml', r'\mag_600mK_0T\mag0T600mK-b3-toto_ppp-cull208-FoFc.dat'), 
(0, r'\mag_1200mK_0T\Na2BaMnP2O8-m3-1200mK-0T.xml', r'\mag_1200mK_0T\Na2BaMnP2O8-m3-1200mK-0T.dat'),

(0.73, r'\mag_600mK_1p2T_Hc\Na2BaMnP2O8-600mK-1p2T-Hc.xml', r'\mag_600mK_1p2T_Hc\Na2BaMnP2O8-600mK-1p2T-Hc.dat'), # Skip?, poor refinement
(1.29, r'\mag_600mK_2p2T_Hc\Na2BaMnP2O8-600mK-2p2T-Hc.xml', r'\mag_600mK_2p2T_Hc\Na2BaMnP2O8-600mK-2p2T-Hc.dat'),
(2.17, r'\mag_600mK_3p8T_Hc\Na2BaMnP2O8-600mK-3p8T-Hc.xml', r'\mag_600mK_3p8T_Hc\Na2BaMnP2O8-600mK-3p8T-Hc.dat'),

(1.06, r'\mag_600mK_1p8T_Hab\Na2BaMnP2O8-m3-600mK-1p8T-Hab.xml', r'\mag_600mK_1p8T_Hab\Na2BaMnP2O8-m3-600mK-1p8T-Hab.dat'),
(1.8,  r'\mag_600mK_3T_Hab\Na2BaMnP2O8-m3-600mK-3T-Hab.xml', r'\mag_600mK_3T_Hab\Na2BaMnP2O8-m3-600mK-3T-Hab.dat'),
(2.66, r'\mag_600mK_4p4T_Hab\Na2BaMnP2O8-m3-600mK-4p4T-Hab.xml', r'\mag_600mK_4p4T_Hab\Na2BaMnP2O8-m3-600mK-4p4T-Hab.dat'),
]
# Hab is field along (1 -1 0) direction
FIELD_DEFINITIONS = {'Hab':np.array([np.sqrt(3)/2, -0.5, 0]), 'Hc':np.array([0,0,1])}

# From [Kim et al. 2022, JPhysCondMat 34, 475803]
# Msat = 5.6 μB/Mn from CW fit
# Msat = 4 μB/Mn from magnetizatio nemasurements
# Mc  = [0.1815, 0.3218, 0.5415]*Msat @ B = [1.2, 2.2, 3.8] T
# Mab = [0.2662, 0.4500, 0.6647]*Msat @ B = [1.8, 3.0, 4.4] T
# 
# Mc  = [0.726 , 1.2872, 2.166 ] muB
# Mab = [1.0648, 1.8   , 2.6588] muB



OVF_HEADER = """# OOMMF OVF 2.0
#
# Segment count: 000001
#
# Begin: Segment
# Begin: Header
#
# Title: SPIRIT NBMPO magnetic moments by  Michal
#
# Desc: -
#
# valuedim: 3   ## field dimensionality
# valueunits: none none none
# valuelabels: spin_x spin_y spin_z
#
## Fundamental mesh measurement unit. Treated as a label:
# meshunit: nm
#
# xmin: -0.1
# ymin: 0
# zmin: 0
# xmax: 0.2
# ymax: 0.173205
# zmax: 0.9
#
# meshtype: rectangular
# xbase: 0
# ybase: 0
# zbase: 0
# xstepsize: 0.1
# ystepsize: 0.0866025
# zstepsize: 0.1
# xnodes: 3
# ynodes: 3
# znodes: 10
#
# End: Header
#
# Begin: Data Text"""

OVF_FOOTER = """# End: Data Text
# End: Segment"""

def get_fraction(kz):
    allowed_denominators = [0, 2, 3, 4, 6, 12]
    fract = Fraction(kz).limit_denominator(12)

    if fract.denominator in allowed_denominators:
        ret = fract
    else:
        ret = kz

    return ret

def get_conditions(filepath: str) -> list[float]:
    """Extract conditions of measurements from the filename."""

    entries = filepath.split('\\')[-2].split('_')

    if len(entries) == 4:
        T = int(entries[1][:-2])
        H = float(entries[2][:-1].replace('p','.'))
        Hstr = entries[3]
        Hdir = FIELD_DEFINITIONS[entries[3]]
    elif len(entries) == 3:
        T = int(entries[1][:-2])
        H = float(entries[2][:-1].replace('p','.'))
        Hstr = ''
        Hdir = np.array([0,0,0])
    else:
        raise ValueError(f"Unexpected filename format: {filepath}")
        
    return T, H, Hdir, Hstr

def parse_xml(filename) -> dict:
    # Get main parameters of the structure
    mtp_result = ET.parse(filename).getroot()
    qz = float(mtp_result.find('sample//propagationVector').get('qz'))

    MN = mtp_result.find("sample//atom[@label='MMN2']")
    C123 = [float(MN.find('magneticMoment').get(rr)) for rr in ['Rx', 'Ry', 'Rz']]
    dC123 = [float(mtp_result.find("fitParameters//atom[@label='MMN2']//sigmaMagneticMoment").get(rr)) for rr in ['Rx', 'Ry', 'Rz']]
    ix123 = [float(irrep.get('Ix')) for irrep in MN.find('.//basisVectors')]
    phase = float(MN.find('magneticMoment').get('phase'))

    Cm, C0, Cp = [C123[i] for i in np.argsort(ix123)]
    dCm, dC0, dCp = [dC123[i] for i in np.argsort(ix123)]

    fit = mtp_result.find("fitParameters//results")

    return dict(
        qz=qz,
        C0=C0,
        dC0=dC0,
        Cm=Cm,
        dCm=dCm,
        Cp=Cp,
        dCp=dCp,
        phase=phase,
        fit=fit
)

def mm(Q: np.ndarray, C0: float, Cm: float, Cp:float, r: np.ndarray, phase: float) ->np.ndarray:
    """Determine the magnetic moment of an ion at position `r`,
    given the coefficients `Ci` of irreps and global `phase`.
    
    Parameters
    ----------
    Q: np.ndarray (3,)
        Modulation vector in r.l.u.
    C0, Cm, Cp: float
        Mode amplitudes.
    r: np.ndarray (3,)
        Position of the ion in r.l.u. r=[u,v,w].
    phase: float
        Global phase shift, exp(2 pi i phi), so in fraction of 2 pi units.  
    """
    phi0 = np.array([0,0,1])
    phiP = np.array([0.612+0.354j,0.707j,0])
    phiM = np.conj(phiP)

    rr =  (C0*phi0 + Cm*phiM + Cp*phiP)*np.exp(2*np.pi*1j * (np.dot(Q, r) + phase) )
    return rr.real

def make_table(results) -> str:
    for rr in results:
        filename = PATH_XML + rr[1]
        result = parse_xml(filename)
        T, H, Hdir, Hstr = get_conditions(filename)

        Hs = {
            '': '-',
            'Hc': '$[001]$', 
            'Hab': f'$[1\\bar{{1}}0]$'
        }[Hstr]
        
        print(f" & ".join([
            f"{T}", 
            f"{H}", 
            f'{Hs}', 
            f"{get_fraction(np.abs(result['qz']))}",
            f"{result['C0']:.2f}\,({result['dC0']:.2f})",
            f"{result['Cm']:.2f}\,({result['dCm']:.2f})",
            f"{result['Cp']:.2f}\,({result['dCp']:.2f})",
            f"{result['phase']:.3f}",
            f"{rr[0]}",
            f"{float(result['fit'].get('Rf')):.2f}",
            ""
            ]), end=" \\\\ \n")
        
def make_mxyz_file(filename: str, magnetization: tuple[float,float,float]) -> None:
    """Create a file with the magnetic moments of the ions in the unit cell.
    OVF format uses xyz coordinate system for the magnetic moments."""
    result = parse_xml(PATH_XML+filename)
    T, H, Hdir, Hstr = get_conditions(filename)

    k = np.array([1/3, 1/3, result['qz']])
    C0 = result['C0']
    Cm = result['Cm']
    Cp = result['Cp']
    phase = result['phase']


    ret = [OVF_HEADER]

    # ret = ["# a = [ 5.37290000e+00,  0.00000000e+00,  0.00000000e+00]"]
    # ret.append("# b =[-2.68645000e+00,  4.65306789e+00,  0.00000000e+00]")
    # ret.append("# c =[ 4.34467945e-16,  7.52520555e-16,  7.09540000e+00]")

    # ret.append("#" + "\t".join(["u", "v", "w", "Mx", "My", "Mz"]))
    HH = np.array([
        [1, np.cos(np.pi/3), 0],
        [0, np.sin(np.pi/3), 0],
        [0, 0, 1]
    ])
    moments = []
    for rz in np.arange(0,10):
        for ry in np.arange(0,3):
            for rx in np.arange(0,3):
                r = np.array([rx, ry, rz])
                m = HH @ mm(k, C0, Cm, Cp, r, phase)
                m += np.array(magnetization)*Hdir

                moments.append(m)
                # ret.append("\t".join([str(x) for x in [rx,ry,rz,m[0],m[1],m[2]]]))  # DEBUG
                ret.append("".join([f"{x:>22.12f}" for x in [m[0],m[1],m[2]]]))

    ret.append(OVF_FOOTER)

    moments_norm = np.linalg.norm(moments, axis=1)
    Mav = np.average(moments_norm)
    sig_Mav = np.var(moments_norm)
    print(f'\tAverage moment: {Mav:.3f} +- {sig_Mav:.3f} μB')
                
    return "\n".join(ret)

def make_vesta_file(filename: str, magnetization: tuple[float,float,float], template: str):
    """Create a VESTA file from template dedicated to visualise NBMPO structure"""
    result = parse_xml(filename)
    T, H, Hdir, Hstr = get_conditions(filename)

    k = np.array([1/3, 1/3, result['qz']])
    C0 = result['C0']
    Cm = result['Cm']
    Cp = result['Cp']
    phase = result['phase']

    vectr_template = """   {N}    {vx:.5f}    {vy:.5f}    {vz:.5f} 0
    6   1    {u}    {v}    {w}
 0 0 0 0 0"""
    vectt_template = "   {N}  1.0 {R}   {G}   {B} 1"

    vectr_block = []
    vectt_block = []


    counter = 0
    for rw in np.arange(0,12+1):
        for rv in np.arange(0,2+1):
            for ru in np.arange(0,2+1):
                counter += 1

                r = np.array([ru, rv, rw])
                m = mm(k, C0, Cm, Cp, r, phase)
                m += np.array(magnetization)*Hdir

                R, G, B = 255, 0, 0
                # import colorsys
                # theta = np.arctan2(np.sqrt(m[0]**2+m[1]**2), m[2])
                # hue = (theta/np.pi + 1)
                # RGB  = colorsys.hsv_to_rgb(hue, 1, 1)
                # R, G, B = [int(255*x) for x in RGB]

                vectr_block.append(vectr_template.format(N=counter, vx=m[0], vy=m[1], vz=m[2], u=ru,v=rv,w=rw))
                vectt_block.append(vectt_template.format(N=counter, R=R, G=G, B=B))

    # print('\n'.join(vectt_block))


    return template.format(VECTOR_DEFINITIONS='\n'.join(vectr_block), VECTOR_STYLE_DEFINITIONS='\n'.join(vectt_block))


@dataclass
class EntryHKL:
    h: int
    k: int
    l: int
    F2c: float
    F2o: float
    F2s: float
    stat: str

    def __post_init__(self):
        if self.F2o < 0:
            self.F2o = 0.0

def load_data(filename):
    """
    Load data from a file.
    The file should contain two columns: Fo and Fc.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    ret = []
    for dd in lines:
        vals = dd.split()
        if len(vals) != 6:
            continue 

        if vals[-1][-1] not in ['0','1','2','3','4','5','6','7','8','9']:
            continue
        
        ee = EntryHKL(
            h=int(float(vals[0])),
            k=int(float(vals[1])),
            l=int(float(vals[2])),
            F2c=float(vals[3]),
            F2o=float(vals[4]),
            F2s=float(vals[5]),
            stat='?'
        )

        ret.append(ee)

    return ret

def plot_FoFc_inset(data: list[EntryHKL], R: float, ax: 'Axis' = None, plot_opts: dict = {}) -> 'Figure':
    color = 'tab:blue'

    if ax is None:
        margins = dict(left=0.15, right=0.95, bottom=0.18, top=0.95)
        fig, ax = plt.subplots(figsize=(7/2.52, 7/2.52), gridspec_kw=margins)
    else:
        fig = ax.get_figure()
    
    # Data
    F2o  = np.array([entry.F2o for entry in data])
    F2c  = np.array([entry.F2c for entry in data])
    dF2o = np.array([entry.F2s for entry in data])

    # Prettifiers
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.set_xlabel(f'$F_c^2$')
    if plot_opts['ylabel']:
        ax.set_ylabel(f'$F_o^2$')

    min_lim = min(min(F2c), min(F2o)) * 0.95 * 0
    max_lim = max(max(F2c), max(F2o)) * 1.1
    # print(min_lim, max_lim)
    ax.set_xlim(min_lim, max_lim)
    ax.set_ylim(min_lim, max_lim)

    labels = [f'$R_1$={R:.2f}', f'$N_\mathrm{{refl}}=${len(data)}']
    # for n,label in enumerate(labels):
    ax.text(0.02, 0.98, labels[0], transform=ax.transAxes, fontsize=12, zorder=100, ha='left', va='top')
    ax.text(0.98, 0.02, labels[1], transform=ax.transAxes, fontsize=12, zorder=100, ha='right', va='bottom')

    ax.text(-0.02, 1.36, plot_opts['letter'], transform=ax.transAxes, ha='right', va='top',
            fontsize=12)
    ax.text(0, 1.36, plot_opts['label'], transform=ax.transAxes, ha='left', va='top',
            fontsize=12)

    # Plotting
    ax.plot(F2c, F2o, 'o', ms=5, c=color, mec='black', mew=0.5, zorder=50)
    ax.errorbar(F2c, F2o, yerr=dF2o, fmt='.', capsize=2, color=color, zorder=0)

    ax.plot([0, max_lim], [0, max_lim], c='black', lw=0.75)

    # This to show N intervals with same number of reflections.
    # F2o.sort()
    # print(len(F2o))
    # F2o_split = np.split(F2o, 4)
    # for F2o_chunk in F2o_split[:-1]:
    #     FF = max(F2o_chunk)
    #     ax.plot([0, FF], [FF, FF], c='gray', lw=1, ls='--')
    #     ax.plot([FF, FF], [0, FF], c='gray', lw=1, ls='--')


    return fig

def plot_FoFC(results):
    result_filenames = [r[1] for r in results]

    margins = dict(left=0.07, right=0.99, bottom=0.12, top=0.9, wspace=0.3, hspace=0.85)
    mosaic = [
        ['0', '2', '3', '4'],
        ['1', '5', '6', '7'],
    ]
    plot_opts = [
        {'ylabel': True,  'letter': 'a)', 'label': '$T$=600 mK,\n$H$=  0'},
        {'ylabel': False, 'letter': 'e)', 'label': '$T$=1200 mK,\n$H$= 0'},
        {'ylabel': True,  'letter': 'b)', 'label': '$T$=600 mK,\n$H$= 1.2 T $\parallel [001]$'},
        {'ylabel': False, 'letter': 'c)', 'label': '$T$=600 mK,\n$H$= 2.2 T $\parallel [001]$'},
        {'ylabel': False, 'letter': 'd)', 'label': '$T$=600 mK,\n$H$= 3.8 T $\parallel [001]$'},
        {'ylabel': False, 'letter': 'f)', 'label': f'$T$=600 mK,\n$H$= 1.8 T $\parallel [1\\bar{{1}}0]$'},
        {'ylabel': False, 'letter': 'g)', 'label': '$T$=600 mK,\n$H$= 3 T $\parallel [1\\bar{{1}}0]$'},
        {'ylabel': False, 'letter': 'h)', 'label': '$T$=600 mK,\n$H$= 4.4 T $\parallel [1\\bar{{1}}0]$'},
    ]

    # for n, letter in zip([item for row in mosaic for item in row], 'abcdefgh'):
    #     n = int(n)
    #     T, H, Hdir, Hstr = get_conditions(result_filenames[int(n)])
    #     field_value_ss = f'{H:.1f}' if H != 0 else '0'
    #     field_dir_ss = {'':'', 'Hc':'$\parallel [001]$', 'Hab':f'$\parallel [1\\bar{{1}}0]$'}[Hstr]
    #     plot_opts[n]['label'] = f'T = {T} mK \nH = {field_value_ss} T {field_dir_ss}'
    #     plot_opts[n]['letter'] = f'{letter})'


    fig, axs = plt.subplot_mosaic(mosaic=mosaic, figsize=(21/2.52, 10/2.52), gridspec_kw=margins)

    for n, res in enumerate(results):
        filename = PATH_XML + res[1]
        result = parse_xml(filename)
    
        plot_FoFc_inset(load_data(PATH_XML + res[2]), float(result['fit'].get('Rf')), ax=axs[f'{n}'], plot_opts=plot_opts[n])
        
        fig_single = plot_FoFc_inset(load_data(PATH_XML + res[2]), float(result['fit'].get('Rf')), ax=None, plot_opts=plot_opts[n])
        fig_single.savefig(PATH_XML + res[2].replace('dat', 'png'), dpi=300)

    return fig

def print_structure(data):
    # filename: str, magnetization: tuple[float,float,float], template: str):
    """Create a VESTA file from template dedicated to visualise NBMPO structure"""
    magnetization = data[0]
    result = parse_xml(PATH_XML+data[1])
    T, H, Hdir, Hstr = get_conditions(data[1])

    k = np.array([1/3, 1/3, result['qz']])
    C0 = result['C0']
    Cm = result['Cm']
    Cp = result['Cp']
    phase = result['phase']

    print('Printing result:', data[1])
    print(f'{k=}')
    print(f'{C0=}')
    print(f'{Cm=}')
    print(f'{Cp=}')
    print(f'{phase=}')

    ret = ''
    counter = 0
    for rw in np.arange(0,12+1):
        ret += f"{rw=} "
        for rv in np.arange(0,0+1):
            ret += f"{rv=} "
            mmm = []
            for ru in np.arange(0,2+1):
                ret += f"{ru=} "
                counter += 1

                r = np.array([ru, rv, rw])
                m = mm(k, C0, Cm, Cp, r, phase)
                m += np.array(magnetization)*Hdir
                mmm.append(m)
                ret += str(m) + ' '

            angs = [np.degrees(np.arctan2(m[2], np.sqrt(m[0]**2+m[1]**2))) for m in mmm]
            print(angs, max(angs)-min(angs))

        ret += '\n'

    print(ret)

if __name__ == "__main__":
    ### Print structure for debugging
    # print_structure(results[0])

    ### Make Fo-Fc plot for all results
    fig = plot_FoFC(results)
    fig.savefig(PATH_XML + r'\NBMPO-FoFc-all.png', dpi=300)
    fig.savefig(PATH_XML + r'\NBMPO-FoFc-all.pdf')

    ### Export to OVF file
    # for magnetization, filename, _ in results:
    #     filename_out = filename.split('\\')[-1]
    #     filename_out = filename_out.replace('.xml', '.ovf')
    #     filename_out = PATH_OVF + '\\' + filename_out
    #     # print(filename)
    #     # print(PATH_OVF+filename_out)
    #     print(f'Making file: {filename_out}')

    #     mxyz = make_mxyz_file(filename, magnetization)
    #     with open(filename_out, 'w') as f:
    #         f.write(mxyz)

    ### Export to VESTA file
    # vesta_template = None
    # with open(r'C:\GDrive\PostDoc-Juelich\7_NBMP-Nikolaos\6_D23-finalized\NBMPO_template.vesta', 'r') as ff:
    #     vesta_template = ff.read()

    # for magnetization, filename, _ in results:
    #     filename_out = filename.split('\\')[-1]
    #     filename_out = filename_out.replace('.xml', '.vesta')
    #     filename_out = PATH_VESTA + '\\' + filename_out
    #     # print(filename)
    #     # print(PATH_OVF+filename_out)
    #     print(f'Making file: {filename_out}')

    #     vesta = make_vesta_file(PATH_XML+filename, magnetization, vesta_template)
    #     with open(filename_out, 'w') as f:
    #         f.write(vesta)