def virial_fit_residuals2(beta_eb_lower, beta_eb_upper, eb, density, potential):
    import subprocess 

    beta_eb = np.zeros(3)
    beta_eb[0] = beta_eb_lower
    beta_eb[2] = beta_eb_upper
    beta_eb[1] = beta_eb[0] + ((beta_eb_upper - beta_eb_lower) / 2)
   
    b2 = np.zeros(3)
    b3 = np.zeros(3)

       # Use virial coeffcients to find T and mu0 with virial coeffcients
    def v_residuals(params, x, b2, b3, eps_data, y=None):
        T = params['T'].value
        mu0 = params['mu0'].value
        #alpha = params['alpha'].value
        #dens_fit = ((2/((2*pi*hbar**2)/((mass_li*kb*T)))) * 
        #           (np.exp((1/(kb*T)) * (mu0 - x))) + 
        #           (2*b2*np.exp(2*(1/(kb*T))*(mu0 -x))) + 
        #           (3*b3*np.exp(3*(1/(kb*T))*(mu0 - x))))

        dens_fit = ((2/((2*pi*hbar**2)/((mass_li*kb*T)))) * 
                   (np.log( 1 + np.exp((1/(kb*T)) * (mu0 - x))) + 
                   (2*b2*np.exp(2*(1/(kb*T))*(mu0 -x))) + 
                   (3*b3*np.exp(3*(1/(kb*T))*(mu0 - x)))))
        if y is None:
            return dens_fit
        return (y - dens_fit) 

    eps_data = 1e-4 
  
    # Repeat bisection procedure for 10 iterations
    for i in range(25):
        for i, beb in enumerate(beta_eb):
           #run virial code to get b2 and b3 coeff
           b2_proc	 = subprocess.Popen(["/home/kristian/dev/c/virial/virial",
                                      str(beta_eb[i]), "1"], stdout=subprocess.PIPE) 
           b3_proc = subprocess.Popen(["/home/kristian/dev/c/virial/virial",
                                       str(beta_eb[i]), "2"], stdout=subprocess.PIPE) 

           b2_out = b2_proc.communicate()[0]
           b3_out = b3_proc.communicate()[0]

           b2[i] = float(b2_out)
           b3[i] = float(b3_out)

        beta_eb_fit = unp.uarray(np.zeros(3),np.zeros(3))
        beta_eb_diff = unp.uarray(np.zeros(3),np.zeros(3))
        # fit for upper, mid and lower beta_eb
        for i, beb in enumerate(beta_eb):
            params = lmfit.Parameters()
            params.add('T', value = 50e-9, min=5e-9, max=70e-9)
            params.add('mu0', value = 1.7e-30, min=1.e-30, max=4e-30)
            #params.add('alpha', value = 1.25, min = 1, max=2)
            fit_output = lmfit.minimize(v_residuals, params, 
                                 args=(potential, b2[i], b3[i], eps_data, density))

            beta_fit = 1/(kb * params['T'].value)
            beta_eb_fit[i] = eb * beta_fit
            beta_eb_diff[i] = np.abs(beta_eb_fit[i] - beb)
            #print(lmfit.fit_report(fit_output))

            #print(beta_eb_diff[i]*100)

        if (beta_eb_diff[0] - beta_eb_diff[1]) < (beta_eb_diff[2] -
                beta_eb_diff[1]):
            beta_eb[2] = beta_eb[1]
            beta_eb[1] = beta_eb[0] + ((beta_eb[2] - beta_eb[0]) / 2)
        else:
            beta_eb[0] = beta_eb[1]
            beta_eb[1] = beta_eb[0] + ((beta_eb[2] - beta_eb[0]) / 2)
   
    beta_eb_max = eb * 1/(kb * (params['T'].value+params['T'].stderr))
    beta_eb_err = np.abs(beta_eb_max - beta_eb[1])
    beta_eb_final = eb * 1/(kb * (params['T'].value))
    fit_output = lmfit.minimize(v_residuals, params, 
                                args=(potential, b2[1], b3[1], eps_data, density))
    
    print('Final betaEb: {0:f} +/- {1:f}'.format(beta_eb_final, beta_eb_err))
    print(lmfit.fit_report(fit_output))
    fit = v_residuals(params, potential, b2[1], b3[1], 1)
    return fit, fit_output#, beta_eb[1]  
