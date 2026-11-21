import numpy as np
import matplotlib.pyplot as plt

def Part_1(delta_t=1800,h_perigee=250,T_target=86164,mu=398600,i_array=np.array([5.2,28.5]),plot=True):
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
    print("Part 1:")
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

    if plot:
        ax = plt.figure().add_subplot(projection='3d')
        phi=np.linspace(0, 2 * np.pi, 100)
        
        #plot earth
        theta = np.linspace(0, np.pi, 10)
        phi_earth=np.linspace(0, 2 * np.pi, 10) #to reduce the number of faces
        x = r_earth * np.outer(np.cos(phi_earth), np.sin(theta))
        y = r_earth * np.outer(np.sin(phi_earth), np.sin(theta))
        z = r_earth * np.outer(np.ones(np.size(phi_earth)), np.cos(theta))
        ax.plot_surface(x, y, z)

        #parametric curve of the GEO
        x0=r_target*np.cos(phi)
        y0=r_target*np.sin(phi)
        z0=np.zeros(100)
        ax.plot(x0, y0, z0, label='GEO')

        #parametric curve of the GTO
        for c in range(len(i_array)):
            x1=a_transfer*np.cos(phi)+a_transfer-r0
            b=(r0*r_target)**0.5 #semi minor axis
            print(b)
            y1=b*np.sin(phi)*np.cos(i_array[c])
            z1=b*np.sin(phi)*np.sin(i_array[c])
            ax.plot(x1, y1, z1, label=f'GTO at {round(i_array[c]*180/np.pi, 2)} degrees')

        # Use the same scale on all three axes so the orbit is not distorted.
        limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
        centers = limits.mean(axis=1)
        radius = 0.5 * np.ptp(limits, axis=1).max()
        ax.set_xlim3d(centers[0] - radius, centers[0] + radius)
        ax.set_ylim3d(centers[1] - radius, centers[1] + radius)
        ax.set_zlim3d(centers[2] - radius, centers[2] + radius)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_zlabel("z (km)")
        plt.legend()
        plt.show()
    return 2*delta_V_rendevous,delta_V_12,i_array

def Part_2(delta_V_01,delta_V_12,i_array,m_payload=2200,Isp_array=np.array([235,315]),epsilon=0.12):
    print("\nPart 2:")
    for indice in range(len(i_array)):
        mass_fraction12 = np.exp(-delta_V_12[indice]*1e3/(Isp_array*9.81)) #separation between tugged mass and untugged mass impulses
        mass_fraction01 = np.exp(-delta_V_01/(Isp_array*9.81))
        mprop12 = m_payload*(1-mass_fraction12)*(1-epsilon)/(mass_fraction12-epsilon)
        mtug= mprop12/(1-epsilon) * 1/mass_fraction01
        print(f"\nFor an inclination of {round(i_array[indice]*180/np.pi, 2)} degrees:")
        for c in range(len(Isp_array)):
            print(f"For an Isp of {Isp_array[c]} s the propellant mass for the Hohmann transfer is: {mprop12[c]} kg and the resulting tug total mass is {mtug[c]}")
    return mtug

delta_V_01,delta_V_12,i_array=Part_1()
Part_2(delta_V_01,delta_V_12,i_array)