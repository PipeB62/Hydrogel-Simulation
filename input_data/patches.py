def main():
    import sys
    import numpy as np

    #dar como argumento sigma y directorio del experimento

    sigma = float(sys.argv[1])
    savedir = f"{sys.argv[2]}/input_data"

    table = (
        "#            eps   sigma   a lambda gamma cos(theta)     A       B     p   q  tol\n"
        "#\n"
        f"Mon Mon Mon 1.000 {sigma} 1.5 0.00 0.00 0.0000000000 14.77811 0.5000 4.0 0.0 0.0\n"
        f"Mon Xl  Xl  1.000 {sigma} 1.5 0.00 0.00 0.0000000000 14.77811 0.5000 4.0 0.0 0.0\n"
        "Xl  Xl  Xl   0.000 0.000   0.0 0.00 0.00 0.0000000000 0.000000 0.0000 0.0 0.0 0.0\n"
        "Mon Mon Xl   0.000 0.000   0.0 0.00 0.00 0.0000000000 0.000000 0.0000 0.0 0.0 0.0\n"
        "Mon Xl  Mon  0.000 0.000   0.0 0.00 0.00 0.0000000000 0.000000 0.0000 0.0 0.0 0.0\n"
        f"Xl  Mon Mon 1.000 {sigma} 1.5 0.00 0.00 0.0000000000 14.77811 0.5000 4.0 0.0 0.0\n"
        "Xl  Mon Xl   0.000 0.000   0.0 0.00 0.00 0.0000000000 0.000000 0.0000 0.0 0.0 0.0\n"
        "Xl  Xl  Mon  0.000 0.000   0.0 0.00 0.00 0.0000000000 0.000000 0.0000 0.0 0.0 0.0"
    )

    with open(f"{savedir}/patches.sw","w") as file:
        file.truncate(0)
        file.writelines(table)

if __name__ == "__main__":
    main()