'''
chain.py

Sequentially run scripts:

select_bandpols.py
generate_pcmt.py


'''
import os, sys

thrs = [0.51,  0.52,  0.53,  0.54,  0.55,  0.56,  0.57,  0.58,
        0.59, 0.61,  0.62,  0.63,  0.64,  0.65,  0.66,  0.67,
        0.68,  0.69, 0.71,  0.72,  0.73,  0.74,  0.75,  0.76,
        0.77,  0.78,  0.79, 0.81,  0.82,  0.83,  0.84,  0.85,
        0.86,  0.87,  0.88,  0.89]

# str_thrs = [th for th in sys.argv[1:]]
str_thrs = [str(th) for th in thrs]
print('str_thrs = ')
print('  '.join(str_thrs)) 

for th in str_thrs:    # ['0.51']
    
    # print('python select_bandpols.py  -s EGHIVY -d /data/geodesy/ ' \
    #       '-o pltEGHIVYm' + th + ' -m ' + th + ' -a')
    # os.system('python select_bandpols.py  -s EGHIVY -d /data/geodesy/ ' \
    #               '-o pltEGHIVYm' + th + ' -m ' + th + ' -a')

    # print('python generate_pcmt.py ' \
    #           'pltEGHIVYm' + th + '/selections_st_EGHIVY_m' + th + \
    #           '.txt pcmt_m' + th)
    # os.system('python generate_pcmt.py ' \
    #               'pltEGHIVYm' + th + '/selections_st_EGHIVY_m' + th + \
    #               '.txt pcmt_m' + th)

    print('python compare_pcmt.py pcmt_m' + th)
    os.system('python compare_pcmt.py pcmt_m' + th)

