%%%gradient noise fractal terrain generation

scale = 2; %generate a larger area
npoints = 300; %adjust level of detail

coefficients = randn(npoints, npoints, 2);
coefficients = normgrads(coefficients, npoints);

X = linspace(1, scale, npoints);
Y = linspace(1, scale, npoints);
C = makeshape(npoints);
%Z = makeheightmap(X, Y, coefficients, npoints);
Z = log(abs(makeheightmap(X, Y, coefficients, npoints)));
%[Xtrans, Ytrans, Ztrans] = sphericalcoords(npoints, scale, X, Y, Z);

handl = surf(X,Y,Z);
    set(handl,'edgecolor','none')
    camlight;
%%
for i = 1:npoints*2
    
    coefficients = randn(npoints, npoints, 2);
    coefficients = normgrads(coefficients, npoints);
    Xfreq = linspace(1, scale*(2^(i/50)), npoints);
    Yfreq = linspace(1, scale*(2^(i/50)), npoints);
    
    %%%Uniformly spaced mountain range
    Zfreq = log(abs(makeheightmap(X+i/10, Y+i/10, coefficients, npoints)));
    
    %%%Noise with shape
    %Zfreq = log(abs(C+makeheightmap(X+rand(1), Y+rand(1), coefficients, npoints)));
    
    %%%Slightly less uniform mountains
    %Zfreq = log(abs(makeheightmap(X+log(i), Y+log(i), coefficients, npoints)));
    
    %%%Essentially just noise
    %Zfreq = log(abs(makeheightmap(X+rand(1), Y+rand(1), coefficients, npoints)));
    
    Z = Z + (2^(-i/50))*Zfreq;
    
    %%%Map to sphere (distance not preserved under the projection)
    %[Xtrans, Ytrans, Ztrans] = sphericalcoords(npoints, scale, Z);
    
    if (mod(i,10) == 0)
        handl = surf(X,Y,Z);
        %handl = surf(Xtrans, Ytrans, Ztrans);
        set(handl,'edgecolor','none')
        camlight;
        drawnow
    end
end

function Z = makeheightmap(X, Y, coefficients, npoints)
    Z = zeros(npoints, npoints);
    for c = 1:npoints
        for r = 1:npoints
            Z(c,r) = gradheightmap(X(c),Y(r), coefficients);
        end
    end
end

function val = interp(a0, a1, w)
    %val = (a1 - a0)*(3 - w*2)*w*w + a0;
    val = (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
end

function grads = normgrads(coefficients, npoints)
    grads = zeros(npoints, npoints, 2);
    for c = 1:npoints
        for r = 1:npoints
            mag = norm([coefficients(r,c,1), coefficients(r,c,2)]);
            grads(r,c,1)=coefficients(r,c,1)/mag;
            grads(r,c,2)=coefficients(r,c,2)/mag;
        end
    end
end

function height = gradheightmap(x, y, coefficients)
    x_floor = floor(x);
    y_floor = floor(y);
    
    vi = [x_floor-x, y_floor-y];
    gradi = [coefficients(x_floor,y_floor, 1), coefficients(x_floor,y_floor, 2)];
    n0 = dot(gradi, vi/sqrt(2));
    vi = [1+x_floor-x, y_floor-y];
    gradi = [coefficients(x_floor+1, y_floor, 1), coefficients(x_floor+1, y_floor, 2)];
    n1 = dot(gradi, vi/sqrt(2));
    ix0 = interp(n0, n1, x-x_floor);
    
    vi = [x_floor-x, 1+y_floor-y];
    gradi = [coefficients(x_floor, y_floor+1, 1), coefficients(x_floor, y_floor+1, 2)];
    n0 = dot(gradi, vi/sqrt(2));
    vi = [1+x_floor-x, 1+y_floor-y];
    gradi = [coefficients(x_floor+1, y_floor+1, 1), coefficients(x_floor+1, y_floor+1, 2)];
    n1 = dot(gradi, vi/sqrt(2));
    ix1 = interp(n0, n1, x-x_floor);
    height = interp(ix0, ix1, y-y_floor);
end

function shape = makeshape(npoints)
    shape = zeros(npoints, npoints);
    for c = 1:npoints
        for r = 1:npoints
            shape(c,r) = (1-(.5-c/npoints)^2-(.5-r/npoints)^2)/log(npoints);
        end
    end
end

function M = lintrans(X, Y, a, b, c, d, npoints)
    M = zeros(npoints, npoints);
    T = [a, b; c, d];
    for c = 1:npoints
        for r = 1:npoints
            M(c,r) = [X(c), Y(r)]*T;
        end
    end
end

function [Xtrans, Ytrans, Ztrans] = sphericalcoords(npoints, scale, Z)

    Xtrans=zeros(npoints,npoints);
    Ytrans=zeros(npoints,npoints);
    Ztrans=zeros(npoints,npoints);
    
    for r = 1:npoints
        for c = 1:npoints
            Xtrans(r,c) = Z(r,c)*cos((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Ytrans(r,c) = Z(r,c)*sin((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Ztrans(r,c) = Z(r,c)*cos((pi*r)/(npoints-1));
        end
    end
end