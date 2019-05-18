# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 13:07:55 2018

@author: mbcxajg2
****************************************
    DATA ANALYSIS: DROP 1 TOP VIEW
****************************************

"""

import numpy as np
import math
import matplotlib.pyplot as plt

#======================================
#           FUNCTION 1           
#======================================
# Returns (contact angle/degree, contact angle error/degree) from inputs: mean radius, mean radius error
def get_theta(mean_radius, mean_radius_error):
    """ 
    float, float -> tuple
    Returns (contact_angle, contact_angle_error) from mean radius
    """
    r = mean_radius * 1e-06                 # convert into SI units
    r_error = mean_radius_error * 1e-06     # convert into SI units
    
    # solve for roots of H 
    h_roots = np.roots([math.pi/6, 0, math.pi*3*(math.pow(r, 2))/6, (-7.6*1e-15)])
    # only want real root of H_roots    
    h = h_roots.real[abs(h_roots.imag)<1e-5] 
    theta = (math.pi/2) - math.atan(((r*r)-(h*h))/(2*r*h))
    
    # calculate h_error by varing r by +- r_error and recalculating h 
    # compute h with r = r + r_error
    h_roots_positve_r_error = np.roots([math.pi/6, 0, math.pi*3*(math.pow((r + r_error), 2))/6, (-7.6*1e-15)]) 
    h_positive_r_error = h_roots.real[abs(h_roots_positve_r_error.imag)<1e-5] 
    # compute h with r = r - r_error
    h_roots_negative_r_error = np.roots([math.pi/6, 0, math.pi*3*(math.pow((r - r_error), 2))/6, (-7.6*1e-15)])
    h_negative_r_error = h_roots.real[abs(h_roots_negative_r_error.imag)<1e-5]
    # take the larget h_error obtained
    if abs(h_positive_r_error - h) > abs(h_negative_r_error - h):
        h_error = abs(h_positive_r_error - h)
    else: 
        h_error = abs(h_negative_r_error - h)
        
    # compute dtheta/dr
    dtheta_dr = -(1/(1+math.pow(((math.pow(r,2)-math.pow(h,2))/(2*r*h)),2)))*((1/(2*h))-(h/(2*r*r)))
    # compute dtheta/dh
    dtheta_dh = -(1/(1+math.pow((((r*r)-(h*h))/(2*r*h)),2)))*(-(r/(2*h*h))-(1/(2*r)))   
    # compute theta_error
    theta_error = math.sqrt( (dtheta_dr*dtheta_dr*r_error*r_error) + (dtheta_dh*dtheta_dh*h_error*h_error) )
    
    return (theta*180/math.pi, theta_error*180/math.pi)
    
#=====================================
#           FUNCTION 2
#=====================================
# Returns (speed/um/s, speed_error/um/s) from inputs: 
#   mean radius, previous mean radius, next mean radius, 
#   time at mean radius, time at previous mean radius, time at next mean radius,
#   mean radius error, previous mean radius error, next mean radius error   
def get_speed(r, prev_r, next_r, t, prev_t, next_t, r_err, prev_r_err, next_r_err):
    """ 
    float * 9 -> tuple
    Calculates contact line speed for every mean radius value except first and last mean radius values
    Returns (speed/, speed_error)
    """
    v1 = (r - prev_r) / (t - prev_t)  # speed at an instance before r
    v2 = (next_r - r) / (next_t - t)  # speed at an instance after r
    v = (v1 + v2) / 2  # speed at r
    v1_err = math.sqrt( ( (1 / (t - prev_t))**2 ) * (r_err**2 + prev_r_err**2) )
    v2_err = math.sqrt( ( (1 / (next_t - t))**2 ) * (next_r_err**2 + r_err**2) )
    v_err = math.sqrt( ((1/4) * (v1_err**2)) + ((1/4) * (v2_err**2)) )
    return (v, v_err)
    

#==========================================
#                LOAD DATA
#==========================================
# rx_data[:, 0] -> time/s | rx_data[:, 1] -> radius/um
r1_data = np.loadtxt('Top_view_drop_1_data_run1.txt')    
r2_data = np.loadtxt('Top_view_drop_1_data_run2.txt')    
r3_data = np.loadtxt('Top_view_drop_1_data_run3.txt')    
data_len = len(r1_data)

#==========================================
#                  TIME
#==========================================
time = r1_data[:, 0]        # time values are same for r1_data to r3_data

#===========================================
#     MEAN RADIUS (UM) AND ERROR
#===========================================
# Calculate mean radius values
mean_r = (r1_data[:, 1] + r2_data[:, 1] + r3_data[:, 1]) / 3

# Array to store mean radius error values
mean_r_err = np.zeros(data_len)      

# Calculate mean radius error values
for i in range(data_len):
    mean_r_err[i] = np.std([r1_data[i][1], r2_data[i][1], r3_data[i][1]]) / math.sqrt(2)

#============================================
#     CONTACT ANGLE (degree) AND ERROR
#============================================
# Arrays to store contact angle and error
theta = np.zeros(data_len)   
theta_err = np.zeros(data_len)   

# Calculate contact angle and error
for i in range(data_len):
    theta[i], theta_err[i] = get_theta(mean_r[i], mean_r_err[i])

#================================================
#     CONTACT LINE SPEED (um/s) AND ERROR
#================================================
# Arrays to store contact line speed and  error
u = np.zeros(data_len)        
u_err = np.zeros(data_len)    

# Calculate contact line speed and error
for i in range(data_len):
    if i == 0:  # no speed for first data point
        continue
    elif i == data_len - 1: # no speed for last data point
        continue
    else: 
        u[i], u_err[i] = get_speed(mean_r[i], mean_r[i-1], mean_r[i+1], time[i], time[i-1], time[i+1], mean_r_err[i], mean_r_err[i-1], mean_r_err[i+1])

# Only want arrays with valid speeds
u = u[1:-1]
u_err = u_err[1:-1]    

#=================================================
#         PLOT: MEAN RADIUS (UM) 
#                   VS
#                 TIME (S)
#=================================================
# Plot
plt.errorbar(time, mean_r, yerr=mean_r_err, marker='o', markersize=2, linestyle='None')
plt.ylabel('Mean Radius / um')
plt.xlabel('Time / s')
plt.grid()
plt.title('Drop 1 Top View\n Mean Radius vs Time')
plt.show()

#==================================================
#           ADJUST CONTACT ANGLE VALUES
#==================================================
# Only want contact angles with valid contact line speeds for graph of 
# contact lines speed against contact angle
theta = theta[1:-1]
theta_err = theta_err[1:-1]

#=================================================
#     PLOT: CONTACT LINE SPEED (UM/S) 
#                   VS
#           CONTACT ANGLE (DEGREE)
#=================================================
# Plot
plt.errorbar(theta, u, xerr=theta_err, yerr=u_err, marker='o', markersize=2, linestyle='None')
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View\n Contact Line Speed vs Contact Angle')
plt.show() 

#=================================================
#          QUADRATIC FIT AND ANALYSIS
#=================================================
# Calculate coefficients and error in coefficients 
(q_coeff, q_covr) = np.polyfit(theta, u, 2, cov=True)

# Calculate u values of quadratic fit
q_u = np.polyval(q_coeff, theta)

# Calculate statistically expected errors for the fit
q_u_expected_err = np.sqrt((1/(1744-3))*(np.sum(np.power(u-q_u, 2)))) # this value is a constant
q_u_expected_err = q_u_expected_err * np.ones(len(q_u)) # create an array of this value

# Plot fitted u against theta
plt.errorbar(theta, q_u, yerr=u_err, marker='o', markersize=3)
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Quadratic Fit]\n Contact Line Speed vs Contact Angle')
plt.show()

# Plot residuals with statistically expected errors
plt.errorbar(theta, q_u-u, yerr=q_u_expected_err, marker='o', markersize=2,  linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Quadratic Fit]\n Residuals with Statistically Expected Errors')
plt.show()

# Plot residuals with actual errors
plt.errorbar(theta, q_u-u, yerr=u_err, marker='o', markersize=2, linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Quadratic Fit]\n Residuals with Actual Errors')
plt.show()

# Calculate chi-square value
q_chisq = np.sum(np.power(((u-q_u)/u_err), 2))
q_reduced_chisq = q_chisq/ (1744-3)
print('[Quadratic Fit] Chi-Sq:', q_chisq)
print('[Quadratic Fit] Reduced Chi-Sq:', q_reduced_chisq)

# Comment on Fit
print('[Quadratic Fit] Comment: Reduced Chi-Sq values within acceptable range. Good fit.')

# Show quadratic equation
print('[Quadratic Fit] Fitted Equation: ax^2 + bx +c')

# Show errors on coefficients
print('[Quadratic Fit] a = %f +- %f' % (q_coeff[0], math.sqrt(q_covr[0][0])))
print('[Quadratic Fit] b = %f +- %f' % (q_coeff[1], math.sqrt(q_covr[1][1])))
print('[Quadratic Fit] c = %f +- %f' % (q_coeff[2], math.sqrt(q_covr[2][2])))

#=================================================
#           CUBIC FIT AND ANALYSIS
#=================================================
# Calculate coefficients and error in coefficients 
(c_coeff, c_covr) = np.polyfit(theta, u, 3, cov=True)

# Calculate u values of cubic fit
c_u = np.polyval(c_coeff, theta)

# Calculate statistically expected errors for the fit
c_u_expected_err = np.sqrt((1/(1744-4))*(np.sum(np.power(u-c_u, 2)))) # this value is a constant
c_u_expected_err = c_u_expected_err * np.ones(len(c_u)) # create an array of this value

# Plot fitted u against theta
plt.errorbar(theta, c_u, yerr=u_err, marker='o')
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cubic Fit]\n Contact Line Speed vs Contact Angle')
plt.show()

# Plot residuals with statistically expected errors
plt.errorbar(theta, c_u-u, yerr=c_u_expected_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cubic Fit]\n Residuals with Statistically Expected Errors')
plt.show()

# Plot residuals with actual errors
plt.errorbar(theta, c_u-u, yerr=u_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cubic Fit]\n Residuals with Actual Errors')
plt.show()

# Calculate chi-square value
c_chisq = np.sum(np.power(((u-c_u)/u_err), 2))
c_reduced_chisq = c_chisq/ (1744-4)
print('[Cubic Fit] Chi-Sq:', c_chisq)
print('[Cubic Fit] Reduced Chi-Sq:', c_reduced_chisq) 

# Comment
print('[Cubic Fit] Comment: Not as good as quadratic fit.')

# Show cubic equation
print('[Cubic Fit] Fitted Equation: ax^3 + bx^2 + cx + d')

# Show errors on coefficients
print('[Cubic Fit] a = %f +- %f' % (c_coeff[0], math.sqrt(c_covr[0][0])))
print('[Cubic Fit] b = %f +- %f' % (c_coeff[1], math.sqrt(c_covr[1][1])))
print('[Cubic Fit] c = %f +- %f' % (c_coeff[2], math.sqrt(c_covr[2][2])))
print('[Cubic Fit] d = %f +- %f' % (c_coeff[3], math.sqrt(c_covr[3][3])))

#=================================================
#           LINEAR FIT AND ANALYSIS
#=================================================
# Calculate coefficients and error in coefficients 
(l_coeff, l_covr) = np.polyfit(theta, u, 1, cov=True)

# Calculate u values of linear fit
l_u = np.polyval(l_coeff, theta)

# Calculate statistically expected errors for the fit
l_u_expected_err = np.sqrt((1/(1744-2))*(np.sum(np.power(u-l_u, 2)))) # this value is a constant
l_u_expected_err = l_u_expected_err * np.ones(len(l_u)) # create an array of this value

# Plot fitted u against theta
plt.errorbar(theta, l_u, yerr=u_err, marker='o')
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Linear Fit]\n Contact Line Speed vs Contact Angle')
plt.show()

# Plot residuals with statistically expected errors
plt.errorbar(theta, l_u-u, yerr=l_u_expected_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Linear Fit]\n Residuals with Statistically Expected Errors')
plt.show()

# Plot residuals with actual errors
plt.errorbar(theta, l_u-u, yerr=u_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Linear Fit]\n Residuals with Actual Errors')
plt.show()

# Calculate chi-square value
l_chisq = np.sum(np.power(((u-l_u)/u_err), 2))
l_reduced_chisq = l_chisq/ (1744-2)
print('[Linear Fit] Chi-Sq:', l_chisq)
print('[Linear Fit] Reduced Chi-Sq:', l_reduced_chisq) 

# Comment
print('[Linear Fit] Comment: Not as good as quadratic fit.')

print('[Linear Fit] Fitted Equation: ax + b')

# Show errors on coefficients
print('[Linear Fit] a = %f +- %f' % (l_coeff[0], math.sqrt(l_covr[0][0])))
print('[Linear Fit] b = %f +- %f' % (l_coeff[1], math.sqrt(l_covr[1][1])))

#===================================================
#           DE GENNES FIT AND ANALYSIS
#===================================================
# Calculate coefficients and error in coefficients 
(g_coeff, g_covr) = np.polyfit(theta*theta, u, 1, cov=True)

# Calculate u values of linear fit
g_u = np.polyval(g_coeff, theta*theta)

# Calculate statistically expected errors for the fit
g_u_expected_err = np.sqrt((1/(1744-3))*(np.sum(np.power(u-g_u, 2)))) # this value is a constant
g_u_expected_err = g_u_expected_err * np.ones(len(g_u)) # create an array of this value

# Plot fitted u against theta*theta
plt.errorbar(theta, g_u, yerr=u_err, marker='o')
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [de Gennes Fit]\n Contact Line Speed vs Contact Angle')
plt.show()

# Plot residuals with statistically expected errors
plt.errorbar(theta, g_u-u, yerr=g_u_expected_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [de Gennes Fit]\n Residuals with Statistically Expected Errors')
plt.show()

# Plot residuals with actual errors
plt.errorbar(theta, g_u-u, yerr=u_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [de Gennes Fit]\n Residuals with Actual Errors')
plt.show()

# Calculate chi-square value
g_chisq = np.sum(np.power(((u-g_u)/u_err), 2))
g_reduced_chisq = g_chisq/ (1744-3)
print('[de Gennes Fit] Chi-Sq:', g_chisq)
print('[de Gennes Fit] Reduced Chi-Sq:', g_reduced_chisq) 

# Comment
print('[de Gennes Fit] Comment: Ok Fit.')

print('[de Gennes Fit] Fitted Equation: ax + b')

# Show errors on coefficients
print('[de Gennes Fit] a = %f +- %f' % (g_coeff[0], math.sqrt(g_covr[0][0])))
print('[de Gennes Fit] b = %f +- %f' % (g_coeff[1], math.sqrt(g_covr[1][1])))


#===================================================
#           COX VOINOV FIT AND ANALYSIS
#===================================================
# Calculate coefficients and error in coefficients 
(v_coeff, v_covr) = np.polyfit(theta*theta*theta, u, 1, cov=True)

# Calculate u values of linear fit
v_u = np.polyval(v_coeff, theta*theta*theta)

# Calculate statistically expected errors for the fit
v_u_expected_err = np.sqrt((1/(1744-4))*(np.sum(np.power(u-v_u, 2)))) # this value is a constant
v_u_expected_err = v_u_expected_err * np.ones(len(v_u)) # create an array of this value

# Plot fitted u against theta
plt.errorbar(theta, v_u, yerr=u_err, marker='o')
plt.ylabel('Contact Line Speed / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cox Voinov Fit]\n Contact Line Speed vs Contact Angle')
plt.show()

# Plot residuals with statistically expected errors
plt.errorbar(theta, v_u-u, yerr=v_u_expected_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cox Voinov Fit]\n Residuals with Statistically Expected Errors')
plt.show()

# Plot residuals with actual errors
plt.errorbar(theta, v_u-u, yerr=u_err, marker='o', linestyle='None')
plt.ylabel('Residuals / um/s')
plt.xlabel('Contact Angle / degree')
plt.grid()
plt.title('Drop 1 Top View [Cox Voinov Fit]\n Residuals with Actual Errors')
plt.show()

# Calculate chi-square value
v_chisq = np.sum(np.power(((u-v_u)/u_err), 2))
v_reduced_chisq = v_chisq/ (1744-4)
print('[Cox Voinov Fit] Chi-Sq:', v_chisq)
print('[Cox Voinov Fit] Reduced Chi-Sq:', v_reduced_chisq) 

# Comment
print('[Cox Voinov Fit] Comment: Ok Fit.')

print('[Cox Voinov Fit] Fitted Equation: ax^3 + b')

# Show errors on coefficients
print('[Cox Voinov Fit] a = %f +- %f' % (v_coeff[0], math.sqrt(v_covr[0][0])))
print('[Cox Voinov Fit] b = %f +- %f' % (v_coeff[1], math.sqrt(v_covr[1][1])))


#=================================================
#       COMPARE ALL CHI-SQ VALUES
#=================================================
print('\n[Quadratic Fit] Reduced Chi-Sq:', q_reduced_chisq)
print('[Cubic Fit] Reduced Chi-Sq:', c_reduced_chisq) 
print('[Linear Fit] Reduced Chi-Sq:', l_reduced_chisq) 
print('[de Gennes Fit] Reduced Chi-Sq:', g_reduced_chisq) 
print('[Cox Voinov Fit] Reduced Chi-Sq:', v_reduced_chisq) 




















