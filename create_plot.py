import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os


figheight = 8
Lx = 64
Ly = 64
uscale = 400  # scale for velocity arrows (higher values -> shorter arrows)
n_rods = min(max(Lx, Ly), 32)
S0 = 0.5 * 1.0

# root = os.getcwd()
root = '/Users/tomi/Desktop/ActiveNematics/my_test_name'
directories = ['u', 'E', 'n_dor', 'ss', 'omega', 'p']
filenames = ['u_0000094001.txt', 'E_0000094001.txt', 'n_dor_0000094001.txt',
             'ss_0000094001.txt', 'omega_0000094001.txt', 'p_0000094001.txt']


def read_data_from_file(directory, filename):
    x, y = [], []
    with open(os.path.join(directory, filename), 'r') as file:
        for line in file:
            x_, y_ = [float(i) for i in line.strip().split()]
            x.append(x_,)
            y.append(y_)
    return x, y

def read_single_column_data_from_file(directory, filename):
    data = []
    with open(os.path.join(directory, filename), 'r') as file:
        for line in file:
            x_ = float(line.strip())
            data.append(x_)
    return data

def read_triple_column_data_from_file(directory, filename):
    nx, ny, dor = [], [], []
    with open(os.path.join(directory, filename), 'r') as file:
        for line in file:
            x_, y_, d_ = map(float, line.strip().split())
            nx.append(x_)
            ny.append(y_)
            dor.append(d_)
    return nx, ny, dor,



def ss_plot_function(ss, vmax):
    """Threshold on absolute vlaues of saddle-splay to show defects"""
    return (np.ma.masked_where(np.abs(ss) < vmax, ss)).T
    

def create_plot(ux, uy, p, Ex, Ey, ω, nx, ny, dor, ss):
    n_quiver_scale = (2/3)* n_rods #highervalue -> SHORTER director lines
    ss_scale = 0.1 * S0
    

    X, Y = np.meshgrid(np.arange(Lx),np.arange(Ly))
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(2.25 * figheight * Lx / Ly, figheight))
    extent = [0, (Lx-1), (Ly-1), 0]
    plots_dict = {
            "X" : X,
            "Y" : Y,
            "ax1" : ax1,
            "ax2" : ax2,
            "ax3" : ax3,
            "strain_tensor" : (Ex, Ey),
            "vorticity" : ω}
        
    fig.tight_layout(pad=5.0)
    ax1.set_box_aspect(Ly/Lx)
    ax2.set_box_aspect(Ly/Lx)
    ax3.set_box_aspect(Ly/Lx)
    ax1.set_title('nematic director field')
    ax2.set_title('velocity field')
    ax3.set_title('pressure field')
        
    # quiver plot for velocity field
    ωscale = 0.1 * uscale / Lx
    plots_dict['ures'] = ures = max(1,int(Lx/20)) # how many velocity arrows    
    plots_dict['u_quiv'] = u_quiv = ax2.quiver(
            X.T[::ures,::ures], Y.T[::ures,::ures],
            ux[::ures, ::ures], uy[::ures, ::ures],
            ω[::ures,::ures], # color velocity arrows by vorticity 
            scale=uscale, 
            angles='xy', 
            cmap='cool', 
            clim=[-ωscale,ωscale]
        )
    ax2.set_facecolor('black')
    ax2.set_xlim(*tuple(extent[:2]))
    ax2.set_ylim(*tuple(extent[2:]))
        
    plt.colorbar(
            u_quiv, 
            ax=ax2, 
            fraction=0.046, 
            pad=0.04, 
            label='ω', 
            location='right',
            extend='both'
        )
        
        # quiver plot for director field    
    plots_dict['nres'] = nres = max(1,int(Lx/n_rods)) # how many director lines
    plots_dict['n_quiv'] = ax1.quiver(
            X.T[::nres,::nres], Y.T[::nres,::nres],
            nx[::nres,::nres], ny[::nres,::nres],
            linewidth=.2, 
            headwidth=0, 
            scale=n_quiver_scale,
            pivot='middle', 
            angles='xy'
        )
        
        # degree of order plot
    plots_dict['do_plot'] = do_plot = ax1.imshow(
            dor.T / S0,
            vmin = 0,
            vmax = 1,
            cmap = 'Greens_r',
            interpolation='none', 
            extent=extent
        )
    plt.colorbar(
            do_plot, 
            ax=ax1, 
            fraction=0.046, 
            pad=0.04, 
            label=r'$S / S_0$', 
            location='right',
            extend='max'
        )
        
        # saddle-splay markers for defects
    plots_dict['ss'] = ss   
    vmax = ss_scale    
    plots_dict['ss_plot'] = ax1.imshow(
            ss_plot_function(ss, vmax),
            vmin=-vmax, 
            vmax=vmax, 
            cmap='bwr_r',
            interpolation='none', 
            extent=extent
        )
        
        # pressure plot
    plots_dict['p_plot'] = p_plot = ax3.imshow(
            p.T,
            norm=colors.Normalize(vmin=p.min(),vmax=p.max()),
            cmap = 'inferno',
            interpolation='none', 
            extent=extent
        )
    plots_dict['p_cbar'] = p_cbar = plt.colorbar(
            p_plot, 
            ax=ax3, 
            fraction=0.046, 
            pad=0.04, 
            label=r'$p$', 
            location='right',
            extend='max'
        )
        
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
        
    plt.show()
    fig.savefig('image.png', dpi=200)


if __name__ == '__main__':
    
    u_dir, E_dir, n_dor_dir, ss_dir, omega_dir, p_dir = [os.path.join(root, d) 
                                                         for d in directories]

    ux, uy = read_data_from_file(u_dir, filenames[0])
    Ex, Ey = read_data_from_file(E_dir, filenames[1])
    nx, ny, dor = read_triple_column_data_from_file(n_dor_dir, filenames[2])
    ss = np.reshape(read_single_column_data_from_file(ss_dir, filenames[3]), (Lx, Ly))
    ω = np.reshape(read_single_column_data_from_file(omega_dir, filenames[4]), (Lx, Ly))
    p = np.reshape(read_single_column_data_from_file(p_dir, filenames[5]), (Lx, Ly))
    
    ux = np.reshape(ux, (Lx,Ly))
    uy = np.reshape(uy, (Lx,Ly))
    nx = np.reshape(nx, (Lx,Ly))
    ny = np.reshape(ny, (Lx,Ly))
    dor = np.reshape(dor, (Lx,Ly))
    
    
    create_plot(ux, uy, p, Ex, Ey, ω, nx, ny, dor, ss)
    
    
    
