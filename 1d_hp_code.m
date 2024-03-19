%  ------> INPUTS <------

%  Input of number of elements
NELEM = 3;

%  Input of length of the bar
XLEN=100;

%  Read the approximation order
IPVAL=3;

%  Read the linear variation of AE
%  Read the coefficients of Linear EQ.
A0=1;
a1=1;

%  Read the Quadratic variation of distributed spring
C0=1;
C1=1;
C2=1;

% Read the quadratic variation of distributed traction

T0=1;
T1=1;
T2=1;

% Boundary Condition at A [Left-most end]

% Displacement Boundary Condition
% Ua='1';

% Force Boundary Condition at A [Left-most end]
% Fa='1';

% Mixed Boundary Condition at A [Left-Most end]
% Ka='1';
% Del_a='1';

% Read the type of Approximating Order
% Type='1'    ----->LAGRANGE
% Type='2'    ----->HIERARCHIC

% Read the point force at the MID-LENGTH
% F_mid='0'   ------>No Point force
% F_mid='1'   ------>Yes, and value is taken


%  -------> Pre-Processing  <----------

NDOF=IPVAL+1;     %Number of nodes in a element

% Mesh Generation
X=zeros(1,NELEM+1);
h_I=XLEN/NELEM;   %Uniform Mesh
for INOD = 1:NELEM+1
    X(INOD)=(INOD-1)*h_I;
end

% Load Integration Points
NINT=6;
x=[];
x(1)=-0.9324695142031520278123016;
x(2)=-0.6612093864662645136613996;
x(3)=-0.2386191860831969086305017;
x(4)=0.2386191860831969086305017;
x(5)=0.6612093864662645136613996;
x(6)=0.9324695142031520278123016;

%Load Corresponding Weight Functions
w=[];
w(1)=-0.1713244923791703450402961;
w(2)=-0.3607615730481386075698335;
w(3)=-0.4679139345726910473898703;
w(4)=0.4679139345726910473898703;
w(5)=0.3607615730481386075698335;
w(6)=0.1713244923791703450402961;

%Load Shape Functions
%------>   LAGRANGE  <---------
for i=1:NINT
    if IPVAL==1
        N=zeros(2,IPVAL+1);
        N(1,i)=(1-x(i))/2;
        N(2,i)=(1+x(i))/2;
    elseif IPVAL==3
        N=zeros(3,IPVAL+1);
        N(1,i)= x(i)*(x(i)-1)/2;
        N(2,i)= x(i)*(x(i)+1)/2;
        N(3,i)= 1-(x(i)*x(i));
    elseif IPVAL==2
        N=zeros(4,IPVAL+1);
        N(1,i)=((1-x(i))*(3*x(i)-1)*(3*x(i)+1))/16;
        N(2,i)=((1+x(i))*(3*x(i)-1)*(3*x(i)+1))/16;
        N(3,i)=9*((x(i)-1)*(3*x(i)-1)*(x(i)+1))/16;
        N(4,i)=9*((1-x(i))x*(1+x(i))*(3*x(i)+1))/16;
    end
end














