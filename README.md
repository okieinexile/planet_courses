# planet_courses  

This has files that calculate the orbits of planets.  

planet_orb(planet) is the core of the module.  It takes data from the file planetarydata.csv and plots an abstract orbit of the planet.  By "abstract orbit", I mean an orbit of correct eccentricity and semi-major axis that is set with its perihelion on the postive x-axis and the sun at the origin.  The results are stored in data/<planet>_orbit.csv. 

put_in_position(planet) takes information from planetarydata.csv and moves the orbital data from <planet>_orbit.csv into the actual physical coordinates in the ecliptic coordinate system.  During this process a matrix will be created and used. This matrix will be stored in data/<planet>_matrix.csv.  The final result in put into data/<planet>_ecliptic.csv.

relative_position(planet1,planet2) takes data from data/<planet1>_ecliptic.csv and data/<planet2>_ecliptic.csv and uses them to calculate the position of planet1 with planet2 as the origin.  These are stored in a file data/<planet1>_relativeto_planet2_ecliptic.csv.

celestial_coords(planet) takes the data from data/<planet>_relativeto_earth_ecliptic.csv, transforms it into celestial coordinates in the equatorial system and stores them in a file report/<planet>_table.csv.
