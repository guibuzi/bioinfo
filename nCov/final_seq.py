import matplotlib.pyplot as plt

colormap = {'bat-SL': 'tab:green', 'pangolin/GD': 'tab:orange', 'RaTG13': 'tab:blue', 
            'pangolin/GX': 'tab:red', 'SARS-related': 'tab:purple'
}

result_2019_ncov = {
    'bat-SL': {'x_l': 2222, 'x_u': 2447, 'y': 0},
    'pangolin/GD': {'x_l': 22779, 'x_u': 23263, 'y': 0},
    'RaTG13': {'x_l': [1, 2448, 23264], 'x_u': [2219, 22768, 30384], 'y': [0, 0, 0]}
}

result_ratg13 = {
    'bat-SL': {'x_l':[1, 9382], 'x_u':[8336, 10135], 'y': [0, 0]},
    'pangolin/GD': {'x_l':[8336, 10136, 15574, 22785], 'x_u':[9380, 14598, 21235, 29338], 'y': [0, 0, 0, 0]},
    'pangolin/GX': {'x_l':[21238], 'x_u':[22753], 'y': [0]},
    'SARS-related': {'x_l':[14600], 'x_u':[15573], 'y': [0]}
}



fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(16, 10))
ax = axes[0]
for k, v in result_2019_ncov.items():
    ax.hlines(xmin=v['x_l'], xmax=v['x_u'], y=v['y'], label=k, color=colormap[k], linewidth=4)
ax.set_ylim(-0.5, 0.5)
ax.annotate('2222', xy=(2222, 0), xytext=(2222, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('2447', xy=(2447, 0), xytext=(2447, -0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('22779', xy=(22779, 0), xytext=(22779, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('23263', xy=(23263, 0), xytext=(23263, -0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.legend()
ax.yaxis.set_visible(False)

ax = axes[1]
for k, v in result_ratg13.items():
    ax.hlines(xmin=v['x_l'], xmax=v['x_u'], y=v['y'], label=k, color=colormap[k], linewidth=4)
ax.set_ylim(-0.5, 0.5)
ax.annotate('8336', xy=(8336, 0), xytext=(8336, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('9380', xy=(9380, 0), xytext=(9380, -0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('10135', xy=(10135, 0), xytext=(10135, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('14598', xy=(14598, 0), xytext=(14598, -0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('15573', xy=(15573, 0), xytext=(15573, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('21235', xy=(21235, 0), xytext=(21235, -0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.annotate('22753', xy=(22753, 0), xytext=(22753, 0.4), textcoords='data', arrowprops=dict(arrowstyle="->"))
ax.legend()
ax.yaxis.set_visible(False)
fig.savefig('2019-ncov.jpg')