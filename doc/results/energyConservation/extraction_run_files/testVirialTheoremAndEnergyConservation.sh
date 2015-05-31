# For every timestep, extract from 140 to 150 timesteps variables phi, vx, vy vz
for t in {140..150}
    do
        ./lineextract -v phi -timestep $t -istart 0 0 0 -iend 50 50 50 -m 0 -o phi$t -uda particletest1.uda.002
        ./lineextract -v p.vx -timestep $t -istart 0 0 0 -iend 50 50 50 -m 0 -o vx$t -uda particletest1.uda.002
        ./lineextract -v p.vy -timestep $t -istart 0 0 0 -iend 50 50 50 -m 0 -o vy$t -uda particletest1.uda.002
        ./lineextract -v p.vz -timestep $t -istart 0 0 0 -iend 50 50 50 -m 0 -o vz$t -uda particletest1.uda.002
    done

python testVirialTheoremAndEnergyConservation.py > virialTheoryAndEnergyConservationResults.dat
