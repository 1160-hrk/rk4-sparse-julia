FROM julia:1.10

WORKDIR /workspace


# Install Git
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Preinstall some common Julia packages (optional)
RUN julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Plots")'

# Initialize project environment
RUN julia -e 'using Pkg; Pkg.activate(".")'

CMD [ "julia" ]
