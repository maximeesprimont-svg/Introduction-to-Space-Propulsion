import numpy as np

def delta_V(delta_t,h_perigee,T_traget,mu,i):
    r_earth=6371 #in km
    r0=r_earth+h_perigee # in km (mu must have the same unit)
    r_target=(mu(T_traget/(2*np.pi))**2)**(1/3)
    a_transfer=(r0+r_target)/2
    T_transfer=2*np.pi*(a_transfer**3/mu)**0.5
    delta_teta=2*np.pi*delta_t/T_transfer
    T_int=T_transfer*(1-delta_teta/(2*np.pi))
    a_int=(mu(T_int/(2*np.pi))**2)**(1/3)
    r_int=2*a_int-r0
    V_tranfer_perigee=(mu*(2/r0-1/a_transfer))**0.5
    V_int_perigee=(mu*(2/r0-1/a_int))**0.5
    delta_V_rendevous=V_tranfer_perigee-V_int_perigee
    V_transfer_apogee=(mu*(2/r_target-1/a_transfer))**0.5
    V_target=(mu/r_target)**0.5
    delta_V_target=V_target-V_transfer_apogee
    delta_V_deviation=V_target*(2*(1-np.cos(i)))
    return 2*delta_V_rendevous+delta_V_target+delta_V_deviation

print(delta_V(1800,250,86160,3.986004418e14,5.2*np.pi/180))