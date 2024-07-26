def main():
    import sys
    import numpy as np

    #Se da como argumento la longitud centro-patch y directorio del experimento

    l = float(sys.argv[1])
    savedir = f"{sys.argv[2]}/input_data"

    mon = (
        "#monomer\n\n"

        "3 atoms\n"
        "2 bonds\n"
        "1 angles\n\n"

        "Coords\n\n"
        
        "1 0.0 0.0 0.0\n"
        f"2 {-l} 0.0 0.0\n"
        f"3 {l} 0.0 0.0\n\n"
        
        "Types\n\n"
        
        "1 1\n"
        "2 2\n"
        "3 2\n\n"
        
        "Bonds\n\n"
        
        "1 1 1 2\n"
        "2 1 1 3\n\n"
        
        "Angles\n\n"
        
        "1 1 2 1 3"
        )
    
    xl = (
        "#xl\n\n"
        
        "5 atoms\n"
        "4 bonds\n"
        "6 angles\n\n"
        
        "Coords\n\n"
        
        "1 0.0 0.0 0.0\n"
        f"2 0.0 0.0 {l}\n"
        f"3 {np.sqrt(8/9)*l} 0.0 {-(1/3)*l}\n"
        f"4 {-np.sqrt(2/9)*l} {np.sqrt(2/3)*l} {-(1/3)*l}\n"
        f"5 {-np.sqrt(2/9)*l} {-np.sqrt(2/3)*l} {-(1/3)*l}\n\n"
        
        "Types\n\n"
        
        "1 3\n"
        "2 4\n"
        "3 4\n"
        "4 4\n"
        "5 4\n\n"
        
        "Bonds\n\n"
        
        "1 1 1 2\n"
        "2 1 1 3\n"
        "3 1 1 4\n"
        "4 1 1 5\n\n"
        
        "Angles\n\n"
        
        "1 2 2 1 3\n"
        "2 2 2 1 4\n"
        "3 2 2 1 5\n"
        "4 2 3 1 4\n"
        "5 2 3 1 5\n"
        "6 2 4 1 5"
        )
    
    with open(f"{savedir}/monomer.mol","w") as file:
        file.truncate(0)
        file.writelines(mon)

    with open(f"{savedir}/xlinker.mol","w") as file:
        file.truncate(0)
        file.writelines(xl)

if __name__ == "__main__":
    main()