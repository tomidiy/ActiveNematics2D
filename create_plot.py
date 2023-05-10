import matplolitlib.pyplot as plt
import numpy as np

Lx, Ly = u.shape[:2]
X, Y = np.meshgrid(np.arange(Lx),np.arange(Ly))
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(2.25 * figheight * Lx / Ly, figheight))
E = np.empty_like(u)
ω = np.empty((Lx,Ly))
EΩ_from_u(u, E, ω)
extent = [0, (Lx-1), (Ly-1), 0]
plots_dict = {
        "X" : X,
        "Y" : Y,
        "ax1" : ax1,
        "ax2" : ax2,
        "ax3" : ax3,
        "strain_tensor" : E,
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
        u[::ures,::ures,0], u[::ures,::ures,1],
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
nx, ny, do = n_do_from_Q(Q) # get nematic director and degree of order
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
        do.T / S0,
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
plots_dict['ss'] = ss = np.zeros(Q.shape[:2])
get_saddle_splay(Q, plots_dict['ss'])    
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

if do_save_plot:
    fig.savefig(savefileprefix + f'{stepcount:10d}.png'.replace(' ','0'))


