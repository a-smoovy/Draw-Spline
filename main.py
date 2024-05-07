import numpy as np
import matplotlib.pyplot as plt
import pygame
import sys

radius = 10
window_size = (800, 600)

def graph_function(window, function):
    # Calculate points for the function and plot them on the window
    for x in range(800):  # Adjust range based on window width
        y = function(x)
        pygame.draw.circle(window, (0, 0, 255), (x, 600 - y), 2)  # Invert y for screen coordinates

def graph_cubic_spline(window, x_values, y_values, b, c, d):
    """
    Graph the cubic spline interpolation on the Pygame window.

    Parameters:
    - window: Pygame window to draw on.
    - x_values: List or array of x-values defining the interpolation intervals.
    - y_values: List or array of corresponding y-values for the interpolation intervals.
    - b, c, d: Coefficients of the cubic spline interpolation.
    """

    for i in range(len(x_values) - 1):
        x_start, x_end = x_values[i], x_values[i + 1]
        num_points = 100  # Adjust based on desired smoothness
        x_vals = np.linspace(x_start, x_end, num_points)
        y_vals = [evaluate_spline(x_values, y_values, x, b, c, d) for x in x_vals]

        # Convert x and y coordinates to screen coordinates
        

        # Draw the line segment on the window
        pygame.draw.lines(window, (0, 0, 255), False, (x_vals,y_vals), width = 5)

def not_a_knot_spline(x_values, y_values):
    """
    Compute the not-a-knot cubic spline interpolation coefficients.

    This function calculates the coefficients b, c, and d for the not-a-knot cubic spline
    interpolation based on the given x and y values.

    Parameters:
    - x_values: List or array of x-values defining the interpolation intervals.
    - y_values: List or array of corresponding y-values for the interpolation intervals.

    Returns:
    - b, c, d: Coefficients of the not-a-knot cubic spline interpolation.
      - b: Coefficients for linear terms.
      - c: Coefficients for quadratic terms.
      - d: Coefficients for cubic terms.
    """

    interval_amount = len(x_values) - 1
    x_difference = np.diff(x_values)
    slopes = np.diff(y_values) / x_difference
    second_deriv = np.zeros(interval_amount)
    for i in range(1, interval_amount):
        second_deriv[i] = 3 * (slopes[i] - slopes[i - 1])

    lower_diag = np.zeros(interval_amount + 1)
    upper_diag = np.zeros(interval_amount + 1)
    right_vect = np.zeros(interval_amount + 1)
    lower_diag[0] = 1
    upper_diag[0] = 0
    right_vect[0] = 0

    for i in range(1, interval_amount):
        lower_diag[i] = 2 * (x_values[i + 1] - x_values[i - 1]) - x_difference[i - 1] * upper_diag[i - 1]
        upper_diag[i] = x_difference[i] / lower_diag[i]
        right_vect[i] = (second_deriv[i] - x_difference[i - 1] * right_vect[i - 1]) / lower_diag[i]
    lower_diag[interval_amount] = 1
    right_vect[interval_amount] = 0
    c = np.zeros(interval_amount+1)
    b = np.zeros(interval_amount)
    d = np.zeros(interval_amount)
    for j in range(interval_amount - 1, -1, -1):
        c[j] = right_vect[j] - upper_diag[j] * c[j + 1]
        b[j] = slopes[j] - x_difference[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * x_difference[j])
    return b, c, d

def evaluate_spline(x_values, y_values, x, b, c, d):
    """
    Evaluate the cubic spline interpolation at a given point x.

    Parameters:
    - x_values: List or array of x-values defining the interpolation intervals.
    - y_values: List or array of corresponding y-values for the interpolation intervals.
    - x: The point at which to evaluate the spline.
    - b, c, d: Coefficients of the cubic spline interpolation.

    Returns:
    - y: The interpolated y-value at the point x using the cubic spline.
    """

    interval_amount = len(x_values) - 1
    idx = np.searchsorted(x_values, x)

    if idx == 0:
        idx = 1
    elif idx == len(x_values):
        idx -= 1

    y = y_values[idx - 1] + b[idx - 1] * (x - x_values[idx - 1]) + \
        c[idx - 1] * (x - x_values[idx - 1]) ** 2 + \
        d[idx - 1] * (x - x_values[idx - 1]) ** 3
    return y

def draw_circle_at(coords, window, circles_list):
    # Draw a circle at the given coordinates and add it to the list
    pygame.draw.circle(window, (255, 0, 0), coords, radius)
    circles_list.append(coords)

def draw_all_circles(window, circles_list):
    # Draw all circles in the list on the window
    for coords in circles_list:
        pygame.draw.circle(window, (255, 0, 0), coords, radius)

def get_mouse_click_coords():
    # Initialize Pygame
    pygame.init()

    # Set up the window
    
    window = pygame.display.set_mode(window_size)
    pygame.display.set_caption("Draw Circle on Click")

    # List to store circle coordinates
    circles_list = []

    mouse_x_values = []
    mouse_y_values = []

    

    

    # Game loop
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.MOUSEBUTTONDOWN:
                # Get mouse click coordinates
                mouse_x, mouse_y = pygame.mouse.get_pos()

                mouse_x_values.append(mouse_x)
                mouse_y_values.append(mouse_y)

                print(f"Mouse clicked at ({mouse_x}, {mouse_y})")
                draw_circle_at((mouse_x, mouse_y), window, circles_list)

                b, c, d = not_a_knot_spline(mouse_x_values, mouse_y_values)
                print(b,c,d)

                # Graph the cubic spline interpolation
                graph_cubic_spline(window, mouse_x_values, mouse_y_values, b, c, d)

        # Clear the window
        window.fill((255, 255, 255))

        '''graph_function(window, function(x_values=mouse_x_values, y_values=mouse_y_values))'''

        

        # Draw all circles
        draw_all_circles(window, circles_list)

        # Update the display
        pygame.display.update()

    # Quit Pygame
    pygame.quit()
    sys.exit()

# Example usage
get_mouse_click_coords()