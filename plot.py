import numpy
from matplotlib import pyplot as plt


def plot_arrays(time, arr1, arr2, col1, col2, xlabel, ylabel, filename):
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(time, arr1, color=col1, linewidth=0.5)
    if arr2 is not None:
        plt.plot(time, arr2, color=col2, linewidth=0.5)
    fig = plt.gcf()
    fig.savefig(filename, dpi=600)
    plt.cla()


def to_interval(x, interval):
    l = interval[1] - interval[0]
    periods = (x - interval[0]) / l
    return x - numpy.floor(periods) * l


parameters = {}
with open("config.ini", "r") as config:
    for s in config:
        if '#' in s:
            s1 = s.split('#')
        else:
            s1 = s
        if '=' in s1:
            parameter = s1.split('=')[0].strip()
            value = s1.split('=')[1].strip().strip('"').strip("'")
            parameters[parameter] = value
n_steps = int(parameters["n_steps"])
m = int(parameters["m"])
n = int(n_steps / m)
t = numpy.zeros(n)
a, a_c = numpy.zeros(n), numpy.zeros(n)
e, e_c = numpy.zeros(n), numpy.zeros(n)
incl, incl_c = numpy.zeros(n), numpy.zeros(n)
Omega, Omega_c = numpy.zeros(n), numpy.zeros(n)
g, g_c = numpy.zeros(n), numpy.zeros(n)
M, M_c = numpy.zeros(n), numpy.zeros(n)
with open("result.dat", "r") as inp:
    k = 0
    for s in inp:
        if '#' in s:
            continue
        s_list = s.split()
        t[k] = float(s_list[0])
        a[k] = float(s_list[1])
        e[k] = float(s_list[2])
        incl[k] = float(s_list[3])
        Omega[k] = float(s_list[4])
        g[k] = float(s_list[5])
        M[k] = float(s_list[6])
        a_c[k] = float(s_list[7])
        e_c[k] = float(s_list[8])
        incl_c[k] = float(s_list[9])
        Omega_c[k] = float(s_list[10])
        g_c[k] = float(s_list[11])
        M_c[k] = float(s_list[12])
        k = k + 1
if "a_units" in parameters.keys():
    a_units = parameters["a_units"]
    a_ylabel = f"$a$, {a_units}"
    da_ylabel = f"$\\Delta a$, {a_units}"
else:
    a_ylabel = "a"
    da_ylabel = "$\\Delta a$"
if "t_units" in parameters.keys():
    t_units = parameters["t_units"]
    t_xlabel = f"$t$, {t_units}"
else:
    t_xlabel = "$t$"
color1 = "black"
color2 = "red"
diff_color = "green"
plot_arrays(t, a, a_c, color1, color2, t_xlabel, a_ylabel, "semimajor_axis.png")
plot_arrays(t, e, e_c, color1, color2, t_xlabel, "eccentricity", "eccentricity.png")
plot_arrays(t, incl, incl_c, color1, color2, t_xlabel, "$i$, degrees", "inclination.png")
plot_arrays(t, Omega, Omega_c, color1, color2, t_xlabel, "$\\Omega$, degrees", "asc_node.png")
plot_arrays(t, g, g_c, color1, color2, t_xlabel, "$g$, degrees", "pericenter_argument.png")
plot_arrays(t, a - a_c, None, diff_color, None, t_xlabel, da_ylabel, "delta_a.png")
plot_arrays(t, e - e_c, None, diff_color, None, t_xlabel, "$\\Delta e$", "delta_e.png")
plot_arrays(t, incl - incl_c, None, diff_color, None, t_xlabel, "$\\Delta i$, degrees", "delta_i.png")
plot_arrays(t, to_interval(Omega - Omega_c, (-180, 180)), None, diff_color, None, t_xlabel, "$\\Delta \\Omega$, degrees", "delta_Omega.png")
plot_arrays(t, to_interval(g - g_c, (-180, 180)), None, diff_color, None, t_xlabel, "$\\Delta g$, degrees", "delta_g.png")
plot_arrays(t, to_interval(M - M_c, (-180, 180)), None, diff_color, None, t_xlabel, "$\\Delta M$, degrees", "delta_M.png")
