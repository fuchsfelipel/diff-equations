"""
This module implements the second, third, and fourth order
Runge-Kutta Numerical Methods.
3 Values are required for the execution of **any** of the RK:
 - dydx:function   --> The function that will be used for RK
 - y0:float        --> The first result for dydx
 - h:float         --> Step size
"""
# This is a programming formality
if __name__ == '__main__':
    # Ask user for data
    input_equation = input("Type the differential equation:\ndy/dx = ")
    y0 = float(input("Type a value for y0: "))
    h = float(input("Type a value for h: "))

    # Transform the input equation into a function dydx that can be used by the computer
    dydx = eval(f'lambda x, y: {input_equation}') 

# Second Order Runge-Kutta method implementation as a function
# The [parameter]:float notations tells that
# the parameter should be a decimal number
def rk_second_order(x0:float, y0:float, x:float, h:float):
    # Determine how many iterations will be used
    # using the desired X in relation to x0 and the step (h)
    n_iter = (int)((x - x0)/h)

    # This is a programming-related formality
    # The parameter y0 cannot be changed
    y = y0 

    # Iterate through the required steps
    # ------------------ Start Iteration -------------------
    for i in range(x0, (n_iter + 1)): # For 1, 2, 3, ..., n_iter
        # Apply the formulas described on the Runge-Kutta section
        m0 = h * dydx(x0, y)
        m1 = h * dydx(x0 + h, y + m0)

        # Update y value
        y = y + (1/2) * (m0+m1)
    # ------------------ Finish Iteration ------------------
    
    # Return (give as an output) the value of the final y
    return y

# Third Order Runge-Kutta method implementation
def rk_third_order(x0:float, y0:float, x:float, h:float):
    n_iter = (int)((x - x0)/h)
    y = y0

    for i in range(1, (n_iter + 1)):
        m0 = h * dydx(x0, y)
        m1 = h * dydx(x0 + h/2, y + m0/2)
        m2 = h * dydx(x + h, y - m0 +2*m1)

        y = y + (1/6) * (m0 + 4*m1 + m2)
    
    return y
    
# Fourth Order Runge-Kutta method implementation
def rk_fourth_order(x0:float, y0:float, x:float, h:float):
    n_iter = (int)((x - x0)/h)
    y = y0

    for i in range(1, (n_iter + 1)):
        m0 = h * dydx(x0, y)
        m1 = h * dydx(x0 + h/2, y + m0/2)
        m2 = h * dydx(x, h/2, y + m1/2)
        m3 = h * dydx(x + h, y + 2*m1)

        y = y + (1/6) * (m0 + 2*m1 + 2*m2 + m3)
    
    return y
