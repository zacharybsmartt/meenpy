from numpy import *
import matplotlib.pyplot as plt


def complete_stress_state(x,y,Txy):
    """
    Given a stress diagram with stress in the x, y and Tau xy, this function outputs the principal stresses, principal plane angles, max shear stress angles, and the average stress
    """
    # Principal stresses and the maximum shear stress
    sigX = x; sigY = y; tauXY = Txy;
    sigma_avg = (sigX + sigY) / 2.0
    R = sqrt(((sigX - sigY) / 2) ** 2 + tauXY ** 2)
    sigma_1 = sigma_avg + R
    sigma_2 = sigma_avg - R
    tau_max = R
    theta_p = 0.5 * 180 * arctan(2 * tauXY / (sigX - sigY)) / pi
    theta_s = 0.5 * 180 * arctan((-sigX + sigY) / (2 * tauXY)) / pi
    if sigX>sigY:
        theta_p1 = theta_p; theta_p2 = theta_p + 90
    elif sigX<sigY:
        theta_p2 = theta_p; theta_p1 = theta_p + 90
    if tauXY>0:
        theta_s1 = theta_s; theta_s2 = theta_s + 90
    elif tauXY<0:
        theta_s2 = theta_s; theta_s1 = theta_s + 90
    print("Principal stresses and the maximum shear stress:")
    print(" sigma_1 =", sigma_1," \n sigma_2 =", sigma_2," \n tau_max =", tau_max)
    print("Principal plane angles:")
    print(" theta_p1 =",theta_p1, "degrees")
    print(" theta_p2 =", theta_p2, "degrees")
    print("Max shear stress angles:")
    print(" theta_s1 =",theta_s1, "degrees")
    print(" theta_s2 =", theta_s2, "degrees")
    print("Average stress: sigma_avg =", (sigX+sigY)/2)
    print("R (MC radius) =", R)


def triple_strain_gauge_to_FoS(a, b, c):
    """
    Given 3 strain gauges on a material with strains a, b, and c this function computes factors of safety as well as all values related to the process calculation, see print statements for output info
    """
    eps_theta_1 = a; eps_theta_2 = b; eps_theta_3 = c
    theta_1 = 0; theta_2 = -45; theta_3 = -90
    E = 200e9; nu=0.3; S_y = 200e6
    a = array([ [1+cos(2*deg2rad(theta_1)),
                1-cos(2*deg2rad(theta_1)),
                2*sin(2*deg2rad(theta_1))],
                [1+cos(2*deg2rad(theta_2)),
                1-cos(2*deg2rad(theta_2)),
                *sin(2*deg2rad(theta_2))],
                [1+cos(2*deg2rad(theta_3)),
                1-cos(2*deg2rad(theta_3)),
                2*sin(2*deg2rad(theta_3))] ])
    b = array([2*eps_theta_1,
            2*eps_theta_2,
            2*eps_theta_3])
    x = linalg.solve(a,b)
    print("Strain in x-frame:")
    print(" eps_x =", round(x[0],2))
    print(" eps_y =", round(x[1],2))
    print(" gamma_xy =", round(2*x[2],2))

    # principal strains and maximum in-plane shear strain
    epsX=x[0]; epsY = x[1]; gammaXY = 2*x[2]
    eps_avg = (epsX + epsY)/2.0
    R = sqrt( ((epsX-epsY)/2)**2 + (gammaXY/2)**2)
    eps_1 = eps_avg + R
    eps_2 = eps_avg - R
    gamma_max = 2*R
    theta_p = 0.5*180*arctan(2*(gammaXY/2)/(epsX-epsY))/pi
    if epsX>epsY:
        theta_p1 = theta_p; theta_p2 = theta_p + 90
    elif epsX<epsY:
        theta_p2 = theta_p; theta_p1 = theta_p + 90
    theta_s1 = theta_p1 - 45; theta_s2 = theta_s1 + 90;
    print("Average strain: eps_avg =", (epsX+epsY)/2)
    print("R (MC radius) =", round(R,2))
    print("Principal strains and the maximum shear strain:")
    print(" eps_1 =", round(eps_1,2)," \n eps_2 =", round(eps_2,2),"\n gamma_max =", round(gamma_max,2))
    print("Principal plane angles:")
    print(" theta_p1 =",round(theta_p1), "degrees")
    print(" theta_p2 =", round(theta_p2), "degrees")
    print("Max shear strain angles:")
    print(" theta_s1 =",round(theta_s1), "degrees")
    print(" theta_s2 =", round(theta_s2), "degrees")
    # principal stresses and maximum shear stress
    eps_1 = eps_1*1e-6; eps_2 = eps_2*1e-6; gamma_max = gamma_max*1e-6
    sigma_2 = E*(eps_1*nu + eps_2)/(1-nu**2)
    sigma_1 = eps_1*E+nu*sigma_2
    G = E/(2*(1+nu))
    tau_max = gamma_max*G
    sigma_vm = sqrt(sigma_1**2 + sigma_2**2 - sigma_1*sigma_2)
    print("Principal stresses an the maximum in-plane shear stress:")
    print(" sigma_1=",round(sigma_1/1e6,2))
    print(" sigma_2=", round(sigma_2/1e6,2))
    print(" tau_max=", round(tau_max/1e6,2))
    sigma_sorted = sort([sigma_1,sigma_2,0])
    sigma_1 = sigma_sorted[2]
    sigma_2 = sigma_sorted[1]
    sigma_3 = sigma_sorted[0]
    print("sigma_1 = ",sigma_sorted[2])
    print("sigma_2 = ",sigma_sorted[1])
    print("sigma_3 = ",sigma_sorted[0])
    print("Failure related:")
    print(" sigma_vm = ", round(sigma_vm/1e6,2), "FOS", S_y/sigma_vm)
    print(" tau_abs-max = ", round(abs((sigma_1-sigma_3)/1e6/2),2), \
    " FOS=", round(S_y/abs(sigma_1-sigma_3),2))


def Mohrs_circle(sigma_x, sigma_y, tau_xy):
    """
    Given sigma_x, sigma_y, tau_xy, produces and plots a Mohrs circle graph
    """    
    try:
        # Convert str to float
        sigma_x = float(sigma_x)
        sigma_y = float(sigma_y)
        tau_xy = float(tau_xy)
        
        # Principal and maximum shear stress values
        p_stress1 = (sigma_x + sigma_y) / 2 +  (((sigma_x - sigma_y) / 2)**2 + tau_xy **2)**0.5
        p_stress2 = (sigma_x + sigma_y) / 2 -  (((sigma_x - sigma_y) / 2)**2 + tau_xy **2)**0.5
        max_shear = (((sigma_x - sigma_y) / 2)**2 + tau_xy **2)**0.5

        print("\\nThe value of Principal stresses are ",p_stress1 , " and ", p_stress2)
        radius = (p_stress1 - p_stress2) / 2
        center_x = (sigma_x + sigma_y) / 2
        center_y = 0

        X = []
        Y = []
        X_new = []

        for i in linspace( 0 , 2 * pi , 150 ):
            x = radius * cos( i ) 
            x_new = x + center_x
            y = radius * sin( i )

            X.append(x)
            X_new.append(x_new)
            Y.append(y)       

        # Plot vertical axis
        if radius > 0:
            y_limit = 1.5*radius
        else:
            y_limit = p_stress1*0.5
        a = linspace(-y_limit, y_limit,100)
        b = a*0
        plt.plot(b, a, color = 'b', linestyle='dashed')

        # Plot horizontal axis
        if p_stress2 < 0:
            x_limit = p_stress2*1.2

        elif radius == 0 :
            x_limit = -abs(p_stress1)*0.5
        else:
            x_limit = -radius*0.5
        a = linspace(x_limit, p_stress1*1.2, 100)
        b = a*0
        plt.plot(a, b, color = 'b', linestyle='dashed')

        # Plot the mohr circle
        plt.plot(X_new, Y, color = 'k')

        # Plot principal stress
        plt.plot((p_stress1), (0), 'o', color = 'b', label = 'P-stress 1')
        plt.plot((p_stress2), (0), 'o', color = 'c', label = 'P-stress 2')

        # Plot the center of Mohr circle
        plt.plot((center_x), (center_y), 'o', color = 'r', label = 'Center')    

        plt.xlabel('Sigma X')
        plt.ylabel('Sigma Y')
        plt.title('Mohr Circle')
        plt.legend()

        # Points on circle
        point1 = [sigma_x,tau_xy ]
        point2 = [sigma_y,-tau_xy]

        # Points on the x-axis
        point3 = [sigma_x, 0]
        point4 = [sigma_y, 0]

        # Define point's coordinates for line
        x12_values = [point1[0], point2[0]]
        y12_values = [point1[1], point2[1]]  
        x31_values = [point3[0], point1[0]]
        y31_values = [point3[1], point1[1]]
        x42_values = [point4[0], point2[0]]
        y42_values = [point4[1], point2[1]]

        # Plot lines within the circle
        plt.plot(x12_values, y12_values, color = 'g')
        plt.plot(x31_values, y31_values, color = 'g')
        plt.plot(x42_values, y42_values, color = 'g')

        plt.plot()
        plt.show()

    except Exception as e:
        print(e)
        print("Sorry, that's not a valid input !")

Mohrs_circle(-20, -40, 10)
