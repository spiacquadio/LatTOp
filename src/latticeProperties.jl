#=
	In this initial implimentation, the cell strut radius will be determined
	base on the porosity. This is because the SIMP method works by penalising 
	the relative density of an element. The porosity of the cell and relative 
	density of the element have a strong corrilation as each define the
	material properties.
=#

using ForwardDiff

#----- Global variables for the lattice structure -----#
# Structure properties
phi = 45.0 * pi / 180.0;	# Angle of the struts
h_cell = 5.0;				# Height of the cell
h=h_cell;
# Material properties
k_pcm = 0.21;		# PCM conductivity
k_cell = 160.0;		# Cell strut material conductivity
E_cell = 70.0e9;	# Cell strut material elastic modulus
rho_pcm = 814.0;	# Density of pcm
rho_cell = 2700.0;	# Denisty of lattice
cp_pcm = 2180.0;	# Specific heat (at constant pressure)
cp_cell = 895.0;	# Specific heat of lattice
L_pcm = 244000.0;	# Latent heat of melting

tm = 273.15 + 29.0;	# Melting temperature
tr = 10.0;			# Half Mushy zone temperature range (T_mushy = [tm-tr : tm+tr])

#----- Cell Property Getters -----#
#=
	Maps the relative density of the elements in the SIMP method
		to a reasonable range of cell volume fraction using a 
		linear interpolation 
=#
function reldense2volfrac(reldense)
	minVF = 0.02;	# Minimum volume fraction
	maxVF = 0.20;	# Maximum volume fraction

	return minVF + (maxVF - minVF) * reldense;
end

#=
	Gets the volume fraction (V_struct/V_total) of the cell
	I.e, this is 1 - porosity
=#
function getVolfrac(r)
	
	theta = atan(tan(phi) / sqrt(2))
	if celltype == "bcc"
		VolFrac = 2*tan(theta^2* (r^2/h_cell^2) * (pi * (1+4/sin(theta)) -16/3* r/h_cell * (3.137/sin(pi-2*theta) + 4.923/sin(pi/2-theta))));
	elseif celltype == "f2ccz"
		VolFrac = tan(phi)^2*r^2/h_cell^2*(pi*(1+4/sin(phi))-16/3*r/h*(2.935/sin(pi-2*phi)+3.667/sin(pi/2-phi)));
	end
	return VolFrac
end;

#= 
	Determines the radius that gives the target porosity
		using the Newton Raphson method
=#
 function determineRadius(targetVolFrac, guess = 1.0)
	# Objective function
	func(r) = getVolfrac(r) - targetVolFrac
	res = 1.0; tol = 1e-8;
	itr = 0; maxItrs = 100;

	# Newton Raphson method used to find actual radius
	while abs(res) > tol && itr < maxItrs
		res = func(guess)
		g = ForwardDiff.derivative(func, guess)

		guess = guess - (res / g)
	end

	return guess
end

#----- Material Property getters -----#
#= 
	Gets the elastic modulus of the cell
	The only input is the cell struct radius as all other properties
		are currently fixed
=#
function getElasticModulus(r)
	E = [];
	if celltype == "bcc"
		prefactor = 2.0 * sqrt(2) * pi * E_cell * r^2 / h_cell^2;
		postfactor = sin(phi)^2 * cos(phi)
		E = prefactor * (1 + 12*r^2/h_cell^2 * (sin(phi)^2 - 1) * tan(phi)^2) * postfactor;
	elseif celltype == "f2ccz"
		prefactor = 2.0*pi*E_cell*r^2/h_cell^2;
		postfactor = sin(phi)^3;
		E = prefactor * (1 + 12*r^2/h_cell^2 * (sin(phi)^2) * tan(phi)^2) * postfactor;
	end
	return E
end;
function getElasticModulus3(r)
	if celltype == "bcc"
	prefactor = 8.0 * pi * E_cell * (r^2 / h_cell^2);
	postfactor = sin(phi)^3 * tan(phi)^2
	E = prefactor * (1 + 12 * r^2/h_cell^2 * cos(phi)^2) * postfactor;
	elseif celltype == "f2ccz"
	prefactor= pi*E_cell*r^2/h_cell^2;
	postfactor = tan(phi)^2;
	E = prefactor * (1+4*sin(phi)^3*(1 + 12*r^2/h_cell^2 * cos(phi)^2)) * postfactor;
	end
	return E
end;

function getElasticityMatrix(r);
	nu = 0.3;
	E12 = getElasticModulus(r);
	E3 = getElasticModulus3(r);
	G13 = E12 / (2 * (1 + nu));

	return (1 / (1 - nu^2)) .* [E12;   nu*E3; 0;;
							    nu*E3; E3;	   0;;
							    0;     0;      G13*(1-nu^2)];
end;

#= 
	Gets the ETC of the cell
	The only input is the cell struct radius as all other properties
		are currently fixed
=#
function getConductivity(r)
	theta = atan(tan(phi) / sqrt(2));
	epsilon = 1 - getVolfrac(r);
	if celltype == "bcc"

		V4 = (2 * (16 / (3 * sin(pi - 2*theta))) * r^3) - (12 * (sqrt(8) - sqrt(6)) * r^3);
		t4 = (V4 * sin(theta) * cos(theta)^2) ^ (1/3)
		Lstr = (h_cell / (2 * sin(theta))) - (t4 * sqrt((2 / cos(theta)^2) + (1 / sin(theta)^2)))

		R1 = (2 * cos(theta)^2) / (4.277 * sin(theta)* t4 * k_cell);
		R2 = Lstr / (pi * r^2 * k_cell)
		R3 = R1;
		Rstr = 2 * (R1 + R2 + R3)
		Rs = Rstr / 4

		k = (epsilon * k_pcm) + (1 / (h_cell * Rs));
	elseif celltype == "f2ccz"
		σ4 = 1.732;
        σ2 = 2.836;        
        V_2Str = 16/3/sin(pi-2*phi)*r^3; # Volumen wenn sich 2 Streben schneiden
        V_4Str = V_2Str*2-12*(8^0.5-6^0.5)*r^3;
        b_2 = (V_2Str*sin(phi)*cos(phi)^2)^(1/3);
        b_4 = (V_4Str*sin(phi)*cos(phi)^2)^(1/3);
        
        L_box = b_4/2*(2/cos(phi)^2+1/sin(phi)^2)^0.5+b_2/2*(2/cos(phi)^2+1/sin(phi)^2)^0.5;
        L_ges = h/(2*sin(phi))-L_box;
        #eps = 1-2*tan(omega)**2/h**3*(4*pi*(2*L_ges)*r**2+2*V_4Str)
        
        R1 = cos(phi)^2/(σ4*sin(phi)*b_4*k_cell);
        R2 = L_ges/(pi*r^2*k_cell);
        R3 = cos(phi)^2/(σ2*sin(phi)*b_2*k_cell);

        R_lat = (2/(2*(R1+R2+R3)))^(-1);
        
        V_ie = pi*r^2*(h/sin(phi)*2+h);
        h_ie = V_ie*tan(phi)/h^2;
        R_pcm = (h/tan(phi)-h_ie)/(h^2/tan(phi)*epsilon*k_pcm)+h_ie/(h^2/tan(phi)*k_cell);
        
        λ_s = 1/(R_lat*h);
        λ_p = 1/(R_pcm*h);
        k = λ_s + λ_p; # [W/(m*K)]
	end
	return k;
end;
function getConductivity3(r)
	theta = atan(tan(phi) / sqrt(2));
	epsilon = 1 - getVolfrac(r);
	if celltype == "bcc"
		V4 = (2 * (16 / (3 * sin(pi - 2*theta))) * r^3) - (12 * (sqrt(8) - sqrt(6)) * r^3);
		t4 = (V4 * sin(theta) * cos(theta)^2) ^ (1/3)
		Lstr = (h_cell / (2 * sin(theta))) - (t4 * sqrt((1 / cos(theta)^2) + (2 / sin(theta)^2)))

		R1 = (2 * sin(theta)^2) / (1.903 * cos(theta)* t4 * k_cell);
		R2 = Lstr / (pi * r^2 * k_cell)
		R3 = R1;
		R_ges = 1/2*(R1+R2+R3);        
        λ_s = 1/(R_ges*h);
        λ_p = epsilon*k_pcm;
        k = λ_s + λ_p; # [W/(m*K)]
	elseif celltype == "f2ccz"
		σ4 = 1.448;
        σ2 = 1.168;
        V_2Str = 16/3/sin(pi-2*phi)*r^3; # Volumen wenn sich 2 Streben schneiden
        V_4Str = V_2Str*2-12*(8^0.5-6^0.5)*r^3;
        b_2 = (V_2Str*sin(phi)^2*cos(phi))^(1/3);
        b_4 = (V_4Str*sin(phi)^2*cos(phi))^(1/3);
        L_box = b_4/2*(1/cos(phi)^2+2/sin(phi)^2)^0.5+b_2/2*(1/cos(phi)^2+2/sin(phi)^2)^0.5;
        L_ges = h/(2*sin(phi))-L_box;
        #eps = 1-2*tan(omega)**2/h**3*(4*pi*(2*L_ges)*r**2+2*V_4Str)
        R1 = 4/2*sin(phi)^2/(σ4*cos(phi)*b_4*k_cell);
        R2 = L_ges/(pi*r^2*k_cell);
        R3 = sin(phi)^2/(σ2*cos(phi)*b_2*k_cell);
        RZ = h/(pi*r^2*k_cell);
        R_lat = (4/(2*(R1+R2+R3))+1/RZ)^(-1);
        λ_s = tan(phi)^2/(R_lat*h);
        λ_p = epsilon*k_pcm;
        k = λ_s + λ_p; # [W/(m*K)]
	end
	return k;
end;

function getConductivityMatrix(r);
	k12 = getConductivity(r);
	k3 = getConductivity3(r);

	return [k12;0;;0;k3];
end

#=
	Gets the density, apparent specific heat and latent heat of an element
=#
function getRhoCpL(reldense)
	VolFrac = reldense2volfrac(reldense);
	rho = VolFrac * rho_cell + (1 - VolFrac) * rho_pcm
	cp = ((VolFrac * rho_cell) / rho * cp_cell) + (((1 - VolFrac) * rho_pcm) / rho * cp_pcm);
	L = (1 - VolFrac) * rho_pcm / rho * L_pcm
	return rho, cp, L;
end;

function getRho(reldense)
	VolFrac = reldense2volfrac(reldense[1]^p);
	rho = VolFrac * rho_cell + (1 - VolFrac) * rho_pcm
	return rho;
end;

function getCp(reldense)
	VolFrac = reldense2volfrac(reldense[1]^p);
	rho = VolFrac * rho_cell + (1 - VolFrac) * rho_pcm
	cp = ((VolFrac * rho_cell) / rho * cp_cell) + (((1 - VolFrac) * rho_pcm) / rho * cp_pcm);
	return cp;
end;

function getL(reldense)
	VolFrac = reldense2volfrac(reldense[1]^p);
	rho = VolFrac * rho_cell + (1 - VolFrac) * rho_pcm
	L = (1 - VolFrac) * rho_pcm / rho * L_pcm
	return L;
end;


using DelimitedFiles
 function exportCells(reldense)
	h = pos[ncon[boundaries["Domain"]][1][3],2] - pos[ncon[boundaries["Domain"]][1][1],2];

	nCells = nElements;
	cellRadii = determineRadius.(reldense2volfrac.(reldense));

	cellPos = zeros(nCells,2);
	 for i = 1:nCells
		cellPos[i,:] = 0.5 * (pos[ncon[boundaries["Domain"]][i][3],1:2] + pos[ncon[boundaries["Domain"]][i][1],1:2]);
	end

	writeData = [h;nCells;cellRadii;cellPos[:,1];cellPos[:,2]];
	writedlm("/Users/spiacqua/Documents/TESI/Topology_Optimisation/res/dehomogenisation/dehomogenize_data.csv", writeData, ",");	
end;