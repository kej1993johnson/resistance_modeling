function [X DX] = trvlsp( N, n) 
% To call this function: [X DX] = trvlsp(10000,10); % a million steps, 20
% cities
% X will return the traveling sequence and DX the travel distance.
% N steps of simulations; n cities to travel (20 max).
% The cities 1 to 20 locations (x,y) are predefined in "city" below.
% Note how the initial traveling sequence is set randomly.

[C city]=distm(n); %distm defined below; C is a symmetric distance matrix
T = 100; % initial Temp
X = randperm(n); %random sequence of number 1 to n, e.g. [11, 3, 1,20,...19]
                 % random sequence of numbers between 1 and n
%Given sequence X what is total distance
DX = distc( n, X, C ); %distc defined below; return the total trvl distance given trvl sequence
% DO MC for all different sequences until we find minimum
for k = 1 : N
    Y = X;
    i = ceil(n*rand); %randomly pick two and swap, generates rand number between 1
                        % and n, ceil rounds up
                       % change the trvl sequence by randomly swap i & j
    j = ceil(n*rand);
    Y([i j]) = X([j i]); % swap the randomly picked numbers
    DY = distc( n, Y, C ); %new total traveling distance
    if rand < exp((DX-DY)./T) %if true, accept the change!
       X = Y; 
       DX = DY;
    end
    T = T*0.99999; % 
    Te(k)=T;
    Score(k)=DX;
end
    %plot T and score over simulation iterations
    subplot (1,2,1), plot(Score, '*');
    hold on;
    plot(Te, '-');
    hold off;
    %plot travel in 2D
    a=num2str(X');
    b=cellstr(a);
    subplot(1,2,2), plot(city(X,1), city(X,2),'ro-'); % plot the cities in the final seq.
    text(city(X,1), city(X,2), b); %& label the cities on the plot
    disp([DX,X]) %print out the traveling distance and sequence
 %end trvlsp
 
 function S = distc(n,x,C) 
%function that return the total traveling distance given the trvl sequence x
%n is the number of cities and C is the matrix storing pair-wise distance.
S = C(x(n),x(1)); %distance between start and end
for i = 1 : n-1, S = S + C(x(i),x(i+1)); end % end distc

function [C city]=distm(n)%store the list of cities, coord. and pairwise dist
C=zeros(n,n);
%each city is defined by (x,y) coordinates
city=  [13     4;
    18     3;
    15     9;
     3    9;
     1    18;
    19     9;
    13    11;
    17    19;
     7     4;
    12    5;
     6    14;
     2    16;
     4    8;
    10     7;
    20    12;
    14    16;
    11    15;
    13     2;
    13    13;
     2    18];
 for i=1:n
     for j= i+1:n
         C(i,j)=(city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2;
         C(j,i)=C(i,j);
     end
 end


