% Function for calculating amount of water vapour at a certain
% height h [km] , given the ground level value m0

function z = wvprof(m0,h)

	scaleh = 2.35;	% scale height 2.35 km
	
	z = m0 * exp(-h / scaleh);
return
