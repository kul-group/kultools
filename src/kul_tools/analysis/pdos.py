class DOSCARReader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.num_atoms = None
        self.num_energy_points = None
        self.energies = []
        self.dos = []

    def read_doscar(self):
        with open(self.file_path, 'r') as file:
            # Read the header information
            line = file.readline().split()
            self.num_atoms = int(line[0])
            self.num_energy_points = int(line[2])

            # Skip the second line
            file.readline()

            # Read the energy values
            for _ in range(self.num_energy_points):
                self.energies.append(float(file.readline().split()[0]))

            # Read the DOS data
            for _ in range(self.num_atoms):
                atom_dos = []
                for _ in range(self.num_energy_points):
                    atom_dos.append(float(file.readline().split()[1]))
                self.dos.append(atom_dos)

    def write_doscar(self, output_file_path):
        with open(output_file_path, 'w') as file:
            # Write the header information
            file.write(f"{self.num_atoms} 1 {self.num_energy_points}\n")
            file.write("Dummy line\n")

            # Write the energy values
            for energy in self.energies:
                file.write(f"{energy} Dummy value\n")

            # Write the DOS data
            for atom_dos in self.dos:
                for dos_value in atom_dos:
                    file.write(f"Dummy value {dos_value}\n")

# Example usage
doscar_reader = DOSCARReader("path/to/DOSCAR")
doscar_reader.read_doscar()

# Access the energy values
energies = doscar_reader.energies

# Access the DOS data
dos = doscar_reader.dos

# Write the DOS data to a new file
doscar_reader.write_doscar("path/to/output/DOSCAR")