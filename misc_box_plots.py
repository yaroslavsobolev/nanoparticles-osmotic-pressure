import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_the_graph(filename, data_for_plot):
    # data_for_plot = data.T
    deltaloc = 0.32
    locations = []
    for i in range(0,8,2):
        locations.extend([i + 1 - deltaloc, i + 1 + deltaloc])
    locations = range(8)
    fig,ax = plt.subplots(1, figsize=(5,5))
    # swp = sns.swarmplot(data=data_for_plot, size=2.5, color='green', zorder=10)
    xs = []
    ys = []
    for i,d in enumerate(data_for_plot):
        xs.extend([locations[i]]*len(d))
        ys.extend(d)
    df = pd.DataFrame(dict(x=xs,
                           y=ys))
    flatui = ['black', 'black', 'orange', 'orange', 'magenta', 'magenta', 'blue', 'blue']
    # sns.set_palette(flatui)
    # sns.palplot(sns.color_palette())

    sns.swarmplot(x="x", y="y", data=df,size=2.5, zorder=10, order=locations, ax=ax, palette=flatui)
    # ax2.artists[2].set_edgecolor('red')
    # for i, x in enumerate(ax2.artists):
    #     if i == 0 or i == 1:
    #         x.set_facecolor('black')
    #     elif i == 2 or i == 3:
    #         x.set_facecolor('orange')
    #     elif i == 4 or i == 5:
    #         x.set_facecolor('magenta')
    #     elif i == 6 or i == 7:
    #         x.set_facecolor('grey')
        # if i % 2 == 0:
        # mybox.set_facecolor('red')
        # mybox.set_edgecolor('black')
        # mybox.set_linewidth(3)


    bp = ax.boxplot(data_for_plot,
                     0, 'rs', 1, whis='range', patch_artist=True,
                     widths=0.6, positions=locations) #list(range(0,8,1))
    # for i, x in enumerate(data_for_plot):
    #     data = x
    #     plt.plot(np.random.uniform(low=i+1-0.2, high=i+1+0.2, size=len(data)), data,
    #          '.', color='green',
    #          alpha=1, markersize=5, zorder=10)
    thealpha = 0.4
    for i,x in enumerate(bp['boxes']):
            x.set_facecolor('white')
            x.set_edgecolor('grey')
            x.set_alpha(0.3)
    for i,x in enumerate(bp['fliers']):
            x.set_color('grey')
            x.set_alpha(thealpha)
    for i,x in enumerate(bp['caps']):
            x.set_color('grey')
            x.set_alpha(thealpha)
    for i,x in enumerate(bp['whiskers']):
            x.set_color('grey')
            x.set_alpha(thealpha)
    # for i,x in enumerate(bp['boxes']):
    #     if i == 0 or i == 1:
    #         x.set_facecolor('black')
    #     elif i == 2 or i == 3:
    #         x.set_color('orange')
    #     elif i == 4 or i == 5:
    #         x.set_facecolor('magenta')
    #     elif i == 6 or i == 7:
    #         x.set_facecolor('grey')
    #     if i % 2 == 0:
    #         x.set_alpha(1)
    #     else:
    #         x.set_alpha(1)
    for median in bp['medians']:
        median.set(color='grey', linewidth=3, linestyle='--', dashes=[8, 2], alpha=0.5)

    plt.ylabel('Galectin3 puncta / cell')
    ax.set_xticklabels(['before', 'after']*4)
    plt.xlabel('')
    ax.set_ylim([-3, 55])
    # plt.tight_layout()
    fig.savefig('data/galectin_data/{0}.png'.format(filename), dpi=600)
    fig.savefig('data/galectin_data/{0}.eps'.format(filename))
    # plt.show()
    # plt.close(fig)
    # del fig
    # del ax
    # del bp

data = np.genfromtxt("data/galectin_data/MCF7-mAG-gal13.txt", delimiter="\t", skip_header=4)
conditions = ['Control', '100:0', '80:20', 'Siramesine']
data_for_plot = []
for d in data.T:
    data_for_plot.append((d[np.logical_not(np.isnan(d))]))
plot_the_graph('b1', data_for_plot)

data = np.genfromtxt("data/galectin_data/U2OS-mchery-gal3.txt", delimiter="\t", skip_header=4)
conditions = ['Control', '100:0', '80:20', 'Siramesine']
data_for_plot = []
for d in data.T:
    data_for_plot.append((d[np.logical_not(np.isnan(d))]))
data_for_plot.append([np.nan])
plot_the_graph('b2', data_for_plot)