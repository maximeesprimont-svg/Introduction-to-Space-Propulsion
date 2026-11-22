# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 16:19:45 2026

@author: lucie
"""

import numpy as np

#Part 1
#parameters
mu=398600
r_earth=6371
h_perigee=250
T_sideral_day=86164
i_1=5.2*np.pi/180
i_2=28.5*np.pi/180
m=2200
delta_t=30*60
r_0=r_earth+h_perigee

#question 1
def question_1():
    # GTO semi-major axis (Hohmann transfer ellipse)
    r_geo=(mu*(T_sideral_day/(2*np.pi))**2)**(1/3)
    a_gto=(r_0+r_geo)/2
    T_gto = 2*np.pi*np.sqrt((a_gto**3)/mu)

    # Angle covered by GTO during the 30 min offset
    d_theta=360*delta_t/T_gto

    # Intermediate phasing orbit period and semi-major axis
    T_int=T_gto*(1-d_theta/360)
    a_int=(mu*(T_int/(2*np.pi))**2)**(1/3)
    r_a_int=2*a_int-r_0

    # Perigee speeds on GTO and on the intermediate orbit
    v_gto_p=np.sqrt(mu*(2/r_0-1/a_gto))
    v_int_p=np.sqrt(mu*(2/r_0-1/a_int))

    # Phasing ΔV: outbound burn and return burn
    dv_go=v_int_p-v_gto_p
    dv_back=-dv_go
    dv_01=abs(dv_go)+abs(dv_back)

    return {
        "r_geo":r_geo, "a_gto": a_gto, "T_gto": T_gto, "d_theta": d_theta,
        "T_int": T_int, "a_int": a_int, "r_a_int": r_a_int,
        "v_gto_p": v_gto_p, "v_int_p": v_int_p,
        "dv_go": dv_go, "dv_back": dv_back,
        "dv_01": dv_01, "t_01": T_int,
    }

#question 2
def question_2(p):
    r_geo=p["r_geo"]
    v_gto_p=p["v_gto_p"]
    v_gto_a=np.sqrt(mu*(2/r_geo-1/p["a_gto"]))
    v_geo=np.sqrt(mu/r_geo)
    T_geo=2*np.pi*np.sqrt((r_geo**3)/mu)
    e_gto=(r_geo-r_0)/(r_geo+r_0)
 
    return {
        "a_gto": p["a_gto"], "T_gto": p["T_gto"], "e_gto": e_gto,
        "r_p": r_0, "r_a": r_geo,
        "v_gto_p": v_gto_p, "v_gto_a": v_gto_a,
        "a_geo": r_geo, "T_geo": T_geo, "v_geo": v_geo,
    }

#question 3
def question_3(orb):
    # ΔV_1-2 (km/s): circularization + inclination change at apogee, in one burn
    v_gto_a=orb["v_gto_a"]
    v_geo=orb["v_geo"]
    dv12_kou=np.sqrt(v_gto_a**2+v_geo**2-2*v_gto_a*v_geo*np.cos(i_1))
    dv12_cape=np.sqrt(v_gto_a**2+v_geo**2-2*v_gto_a*v_geo*np.cos(i_2))
 
    return {
        "dv12_kou": dv12_kou, "dv12_cape": dv12_cape,
    }

#question 4
def question_4(p, dv12):
    # Total ΔV = phasing + circularization/plane change
    dv_tot_kou=p["dv_01"]+dv12["dv12_kou"]
    dv_tot_cape=p["dv_01"]+dv12["dv12_cape"]
 
    return {
        "dv_tot_kou": dv_tot_kou, "dv_tot_cape": dv_tot_cape,
        "difference": dv_tot_cape-dv_tot_kou,
    }

p=question_1()
orb=question_2(p)
dv12=question_3(orb)
tot=question_4(p, dv12)

print("\n==================== PART 1 ====================")
print("Question 1:")
print(f"a_GTO                 : {p['a_gto']:.1f} km")
print(f"T_GTO                 : {p['T_gto']:.1f} s")
print(f"ΔΘ                    : {p['d_theta']:.2f} deg") 
print(f"T_int                 : {p['T_int']:.1f} s")
print(f"a_int                 : {p['a_int']:.1f} km")
print(f"r_a,int               : {p['r_a_int']:.1f} km")
print(f"v_GTO,p               : {p['v_gto_p']:.4f} km/s")
print(f"v_int,p               : {p['v_int_p']:.4f} km/s")
print(f"ΔV_go                 : {p['dv_go']*1000:.1f} m/s")
print(f"ΔV_back               : {p['dv_back']*1000:.1f} m/s")
print(f"ΔV_0-1                : {p['dv_01']*1000:.1f} m/s")
print(f"t_0-1                 : {p['t_01']:.1f} s")
 
print("\nQuestion 2:")
print(f"a_GTO                 : {orb['a_gto']:.1f} km")
print(f"e_GTO                 : {orb['e_gto']:.4f}")
print(f"a_GEO                 : {orb['a_geo']:.1f} km")
print(f"v_GTO,p               : {orb['v_gto_p']:.4f} km/s")
print(f"v_GTO,a               : {orb['v_gto_a']:.4f} km/s")
print(f"v_GEO                 : {orb['v_geo']:.4f} km/s")
print(f"T_GEO                 : {orb['T_geo']:.1f} s")
 
print("\nQuestion 3:")
print(f"ΔV_1-2 Kourou         : {dv12['dv12_kou']:.4f} km/s")
print(f"ΔV_1-2 Cape           : {dv12['dv12_cape']:.4f} km/s")
 
print("\nQuestion 4:")
print(f"ΔV_tot Kourou         : {tot['dv_tot_kou']:.4f} km/s")
print(f"ΔV_tot Cape Canaveral : {tot['dv_tot_cape']:.4f} km/s")
print(f"Difference            : {tot['difference']*1000:.1f} m/s")


#Part 2
#parameters part 2
g0=9.80665
eps=0.12               
isp={"Hydrazine": 235, "MMH/NTO": 315}   # s
 
#total dV of part 1 (km/s)
dv_tot={"Kourou": tot["dv_tot_kou"], "Cape Canaveral": tot["dv_tot_cape"]}
 
#question 2 (part 2)
def question_2_p2():
    res={}
    for site, dv in dv_tot.items():
        dv_ap = dv12["dv12_kou"] if site == "Kourou" else dv12["dv12_cape"]   
        for prop, I in isp.items():
            c = g0*I
            A = np.exp(p["dv_01"]*1000/c)               
            B = np.exp(dv_ap*1000/c)                    
            M1 = m*(B-1)/(1-B*eps)                       
            M0 = M1*A                                    
            res[(site, prop)] = {"c": c, "ratio": A*B, "prop_frac": 1-1/(A*B),
                                 "M0": M0, "M1": M1,
                                 "m_dry": eps*M1,
                                 "m_prop_ph": M0-M1, "m_prop_ap": M1*(1-eps),
                                 "m_prop": M0-eps*M1,
                                 "m_tug": M0, "m0": m+M0}
    return res
 
#question 3 (part 2)
def question_3_p2():
    res = question_2_p2()
    best = min(res, key=lambda k: res[k]["m0"])
    return res, best
 
q2=question_2_p2()
q3, best=question_3_p2()

print("\n==================== PART 2 ====================")
print("\nQuestion 2:")
for (site, prop), r in q2.items():
    print(f"{site} | {prop}")
    print(f"c = g0*Isp            : {r['c']:.1f} m/s")
    print(f"m_i/m_f               : {r['ratio']:.4f}")
    print(f"m_prop/m_i            : {r['prop_frac']:.4f}")
 
print("\nQuestion 3:")
for (site, prop), r in q3.items():
    print(f"{site} | {prop}")
    print(f"m_prop                : {r['m_prop']:.1f} kg")
    print(f"m_dry                 : {r['m_dry']:.1f} kg")
    print(f"m_tug                 : {r['m_tug']:.1f} kg")
    print(f"m0                    : {r['m0']:.1f} kg")
print(f"\nBest configuration    : {best[0]} + {best[1]} (m0 = {q3[best]['m0']:.1f} kg)")
 

#Part 3
#parameters part 3 (MMH/NTO, room temperature)
OF=2.2                
rho_fuel=874          
rho_ox=1440          
burn_fraction=0.01      
 
#winning configuration of part 2
site_win, prop_win=best
isp_win=isp[prop_win]
m_prop_win=q3[best]["m_prop"]
M0_win=q3[best]["M0"]
M1_win=q3[best]["M1"]
dv12_win={"Kourou": dv12["dv12_kou"], "Cape Canaveral": dv12["dv12_cape"]}[site_win]
 
#question 1 (part 3)
def question_1_p3():
    m_fuel=m_prop_win/(1+OF)
    m_ox=m_prop_win*OF/(1+OF)
    v_fuel=m_fuel/rho_fuel
    v_ox=m_ox/rho_ox
    r_fuel=(3*v_fuel/(4*np.pi))**(1/3)   
    r_ox=(3*v_ox/(4*np.pi))**(1/3)
    return {"m_fuel": m_fuel, "m_ox": m_ox, "v_fuel": v_fuel, "v_ox": v_ox,
            "v_tot": v_fuel+v_ox, "r_fuel": r_fuel, "r_ox": r_ox}
 
#question 2 (part 3)
def question_2_p3():
    c=g0*isp_win
    # mass at the start of the apogee burn (after the phasing burn dV_0-1)
    m_start=M1_win+m                      
    m_end=eps*M1_win+m                     
    m_burn=m_start-m_end
    t_burn_max=burn_fraction*p["T_gto"]
    mdot=m_burn/t_burn_max
    thrust=mdot*c
    return {"c": c, "m_start": m_start, "m_burn": m_burn,
            "t_burn_max": t_burn_max, "mdot": mdot, "thrust": thrust}
 
q1_p3=question_1_p3()
q2_p3=question_2_p3()

print("\n==================== PART 3 ====================")
print(f"\nWinning configuration : {site_win} + {prop_win}")
print(f"m0                    : {q3[best]['m0']:.1f} kg")
print(f"m_prop                : {m_prop_win:.1f} kg")
 
print("\nQuestion 1:")
print(f"O/F                   : {OF:.2f}")
print(f"m_fuel (MMH)          : {q1_p3['m_fuel']:.1f} kg")
print(f"m_ox (NTO)            : {q1_p3['m_ox']:.1f} kg")
print(f"V_fuel                : {q1_p3['v_fuel']:.4f} m^3 ({q1_p3['v_fuel']*1000:.1f} L)")
print(f"V_ox                  : {q1_p3['v_ox']:.4f} m^3 ({q1_p3['v_ox']*1000:.1f} L)")
print(f"V_total               : {q1_p3['v_tot']:.4f} m^3")
print(f"R_sphere fuel tank    : {q1_p3['r_fuel']:.3f} m")
print(f"R_sphere ox tank      : {q1_p3['r_ox']:.3f} m")
 
print("\nQuestion 2:")
print(f"T_GTO                 : {p['T_gto']:.1f} s")
print(f"t_burn,max (1 %)      : {q2_p3['t_burn_max']:.1f} s")
print(f"m before apogee burn  : {q2_p3['m_start']:.1f} kg")
print(f"dV_1-2                : {dv12_win*1000:.1f} m/s")
print(f"c = g0*Isp            : {q2_p3['c']:.1f} m/s")
print(f"m_prop apogee burn    : {q2_p3['m_burn']:.1f} kg")
print(f"mdot                  : {q2_p3['mdot']:.4f} kg/s")
print(f"T_min                 : {q2_p3['thrust']:.1f} N")