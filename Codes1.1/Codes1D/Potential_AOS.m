%%%%%%%%%%%%%%%% This function produces the Modified Gravity V_{AOS} potential %%%%%%%%%

function [Vaos] = Potential_AOS(x,e,rH,l)

%x = linspace(-50,50,100);
%m = 1;
%[Rsch, Rsch2m] = RstarToR(x,m);
%rH = 2; the function has been tested for rH = 1, 2.
%e = 0;

%%%%%%%%%%%% We use the tortoise(r*) to radial(r) conversion routine here %%%%%%%%%%%%%%

[Raos, ~] = RstarToR_AOS(x,rH,e);
%l = 0; the function has been tested for l = 0, 2.

Vaos = zeros(size(x,1),size(x,2));
for i = 1:length(Raos)
    a = (Raos(i)/rH)^(2*e);
    b = 1-(rH/Raos(i))^(1+e);
    c = (1+e*(1+ rH/Raos(i)))* l*(l+1)/(((1+e)^3)* Raos(i)^2);
    d = (e+ (rH/Raos(i))^(1+e))/((Raos(i)^2)* (1+e)^2); 
    Vaos(i) = a*b*(c+d);
end

%plot(Raos,V,'-k')
%xlabel('Radial Coordinate')
%ylabel('Potential')
end

