# illustrate how to make and display the Shepp-Logan ellipse digital image phantom
using MIRT: jim, ellipse_im
x = ellipse_im(100)
jim(x, title="Shepp-Logan Phantom")
