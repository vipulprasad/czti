import pandas as pd
pd.options.display.float_format = '{:,.4f}'.format


folder = input("Enter the folder name(eg:20210726_T04_030T01_9000004588_level2_31503): ")
orbit = folder[38:43]
obs = folder[9:12]+'_'+folder[26:30]
data = pd.read_html('https://www.iucaa.in/~astrosat/czti_dqr/{}/index.html'.format(folder))

if len(data) == 10:
    n = 3
else:
    n = 2

print('\nTelemetry error for Observation {} orbit {}'.format(obs, orbit))

for i in range(1,n+1):

        mode_tab = data[i]
        print('\n'+mode_tab.columns[0][0])
        mode_tab.columns = ['Quad', 'bad', 'total', 'discard']
        mode_tab['ratio'] = mode_tab[mode_tab.columns[3]]/mode_tab[mode_tab.columns[2]] * 100
        mode_tab['spc'] = '    '
        #mode_tab['ratio'].round(decimals=4)
        print(mode_tab[[mode_tab.columns[0], 'spc', 'ratio']].to_string(index=False, header=False))
