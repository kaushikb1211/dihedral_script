import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
data = np.random.normal(0,2,2*20000*20).reshape(2,20000,20).clip(-3.14,3.14)
names = np.ones(20).astype(str); names[[3,5,2,0]] ="ALA"
#import matplotlib.pyplot.colors as colors
def plot_dihedral(data:"np.array with dimensions: (phi,psi),(N_frames),(N_residues-1)",
                  resnames:"3 letter strings of sequence", weights:"reweight data", 
                  residue:"residue type to plot", stride=1):
    
    ###
    dic = {residue:{}}
    arr_i = []
    arr_r = []
    idx_= np.where(resnames[1:-1]==residue)[0]
    ###
    
    for idx in idx_:
        fig, ax = plt.subplots(2, 2, figsize=((8,8)), sharex = True, sharey = True)
        hist_i=ax[0,0].hist2d(data[0,::stride,idx+1].flatten(), data[-1,::stride,idx].flatten(), bins=50,range=[[-3.14,3.14],[-3.14,3.14]], norm=colors.LogNorm(),cmap='jet')
        hist_r=ax[0,1].hist2d(data[0,::stride,idx+1].flatten(), data[-1,::stride,idx].flatten(),weights=weights, bins=50,range=[[-3.14,3.14],[-3.14,3.14]], norm=colors.LogNorm(),cmap='jet')
        hist_init=hist_i[0]/np.sum(hist_i[0])
        hist_re=hist_r[0]/np.sum(hist_r[0])
        xedges=hist_i[1]
        yedges=hist_i[2]

        hist_diff=hist_init - hist_re
        ax[1,0].pcolorfast(xedges, yedges, hist_diff.T,cmap='seismic')
        ax[1,1].pcolorfast(xedges, yedges, -hist_diff.T,norm=colors.LogNorm(),cmap='jet')
        
    
        ###
        dic[residue][f"{idx+2}"] = {"initial":hist_i[0], "final":hist_r[1]}
        arr_i.append(hist_i[0])
        arr_r.append(hist_r[0])
        ###
        
        for axes in ax.flat:
            axes.set_xlabel(r"$\phi$", size = 20)
            axes.set_ylabel(r"$\psi$",size = 20)
        fig.suptitle(f"{residue}:{idx+2}", fontsize=20)
   
    #####
    plt.figure()
    arr_r = np.stack(arr_r).mean(0)
    arr_i = np.stack(arr_i).mean(0)
    plt.imshow( arr_r, norm=colors.LogNorm(), interpolation = "bessel")
    plt.figure()
    plt.imshow( arr_i, norm=colors.LogNorm(), interpolation = "bessel")
    dic[residue]["average"] = dict(zip("initial,final".split(","), [arr_r,arr_i]))
    ###
    
    return dic, arr_r
