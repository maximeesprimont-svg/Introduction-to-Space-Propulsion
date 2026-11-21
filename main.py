import numpy as np

def Part_1(delta_t,h_perigee,T_target,mu,i_array):
    r_earth=6371 #in km
    r0=r_earth+h_perigee # in km (mu must have the same unit)
    r_target=(mu*(T_target/(2*np.pi))**2)**(1/3)
    a_transfer=(r0+r_target)/2
    T_transfer=2*np.pi*(a_transfer**3/mu)**0.5
    delta_teta=2*np.pi*delta_t/T_transfer
    T_int=T_transfer*(1-delta_teta/(2*np.pi))
    a_int=(mu*(T_int/(2*np.pi))**2)**(1/3)
    r_int=2*a_int-r0
    V_transfer_perigee=(mu*(2/r0-1/a_transfer))**0.5
    V_int_perigee=(mu*(2/r0-1/a_int))**0.5
    delta_V_rendevous=V_transfer_perigee-V_int_perigee
    V_transfer_apogee=(mu*(2/r_target-1/a_transfer))**0.5
    V_target=(mu/r_target)**0.5
    i_array*=np.pi/180 # to convert to radians
    delta_V_12=(V_transfer_apogee**2+V_target**2-2*V_transfer_apogee*V_target*np.cos(i_array))**0.5

    print(f"Question 1:\ntotal Delta V for rendezvous: {2*delta_V_rendevous} km/s")
    print(f"elapsed time for rendezvous: {T_int} s")

    print(f"\nQuestion 2:\nsemi-major axis of GTO: {a_transfer} km")
    print(f"velocity at perigee of GTO: {V_transfer_perigee} km/s")
    print(f"Velocity at apogee of GTO: {V_transfer_apogee} km/s")
    print(f"semi-major axis of GEO: {a_int} km")
    print(f"Velocity at GEO: {V_target} km/s")

    for c in range(len(i_array)):
        print(f"\nFor an inclination of {round(i_array[c]*180/np.pi, 2)} degrees:")
        print(f"Question 3:\ntotal Delta V for Hohmann transfer: {delta_V_12[c]} km/s")
        print(f"Question 4:\ntotal Delta V for rendezvous and Hohmann transfer: {2*delta_V_rendevous+delta_V_12[c]} km/s")
    return 2*delta_V_rendevous+delta_V_12

Part_1(1800,250,86164,398400.4418,np.array([5.2,28.5]))