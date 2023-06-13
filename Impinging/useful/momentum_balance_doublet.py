import numpy as np
from useful_equation import mdot2velocity

def doublet_sheet_veloity(type, mdot_1, do_1, mdot_2, do_2, rho_1, rho_2, imping_ang):
    #based on momentum balance between 2 impinging jets either identical
    # or unlike jets
    # Assuming that the jets are built symmetrically
    # So the resulting jet is not aligned with the chamber of the jets are unlike


    beta = imping_ang
    orif_ang = beta/2

    if type == 'like':
        v_1 = mdot2velocity(mdot_1, rho_1, do_1)
        Us = v_1 * np.cos(beta/2 * np.pi/180)

        return [v_1, Us]
    elif type == 'unlike':

        v_1 = mdot2velocity(mdot_1, rho_2, do_2)
        v_2 = mdot2velocity(mdot_2, rho_2, do_2)

        
        num_phi = (mdot_1*v_1 - mdot_2*v_2)*np.tan(orif_ang*np.pi/180)
        den_phi = mdot_1*v_1 + mdot_2*v_2
        phi = np.arctan(num_phi/den_phi) * 180/np.pi

        dynamic_term = (mdot_1*v_1 - mdot_2*v_2)/(mdot_1 + mdot_2)

        phi_rad = phi*np.pi/180
        orif_ang_rad = orif_ang*np.pi/180

        Usx = np.cos(orif_ang_rad)/np.cos(phi_rad) * dynamic_term
        Usy = np.sin(orif_ang_rad)/np.sin(phi_rad) * dynamic_term
        Us = np.sqrt(Usx**2 + Usy**2)
        return [v_1, v_2, Us]
    
