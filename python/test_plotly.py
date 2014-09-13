__author__ = 'clipo'

import plotly
py = plotly.plotly(username='clipo', key='fzp2b14ylv')

trace0 = {'x': [1,2,3,4],
  'y': [10,15,13,17],
  'type': 'scatter',
  'mode': 'markers'}

trace1 = {'x': [2,3,4,5],
  'y': [16,5,11,9],
  'type': 'scatter',
  'mode': 'lines'}

trace2 = {'x': [1,2,3,4],
  'y': [12,9,15,12],
  'type': 'scatter',
  'mode': 'lines+markers'}

response = py.plot([trace0, trace1, trace2])
url = response['url']
filename = response['filename']