from subprocess import run
import numpy as np
import my_python_scripts.KKR_QE_plotting as KQP

temps = np.arange(1.025, 1.55, 0.025)
folder = '.'
output_folder = 'output-Field8-0T-Ba_high_res_stat2'
pair_file = '../inputfiles/pair_Mn_J2Ja_small_Jnn0.11_tripleUC_15mat.UCF'

run('mkdir -p '+output_folder, shell=True, cwd=folder)
run('cp input vampire*.mat '+pair_file+' '+output_folder+'/', shell=True, cwd=folder+'')

for temperature in temps:
    KQP.editing_file(filename=folder+'/input', key='sim:temperature=', value=temperature)
    run('mpirun -np 48 vampire-parallel', shell=True, cwd=folder+'')
    run('cp output '+output_folder+'/output_'+f'{temperature:.3f}', shell=True, cwd=folder+'')
    run('cp log '+output_folder+'/', shell=True, cwd=folder+'')
