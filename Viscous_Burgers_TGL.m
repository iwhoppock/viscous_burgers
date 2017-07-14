
%Please Pick your favourite ALPHA (i.e., diffusion coefficient) and RUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.01;%diffusion coefficient%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we must have a square matrix, so there are equal time and space steps
N = 100;

%time
t_max = 2;
dt = t_max / N;

%initialise matrices:
U = zeros(N,N);     %results
D1 = zeros(N,N);    %first derivative matrix

%Obtain TGL points
TGL = flipud(cos(pi * (0:N-1) / (N-1))'); 

%%%%%BEGIN FILLING D1%%%%%
%The following is the result of equation (7.40) in the book

%Fill using MATLAB indexing: START with the  corners 
D1(1,1) = (1 + (2 * (N-1)^2))/6;
D1(N,N) = -(1 + (2 * (N-1)^2))/6;
D1(1,N) = 0.5 * (-1)^(N-1);
D1(N,1) = -0.5 * (-1)^(N-1);

%second row:
for k = 2:(N-1)
    D1(1,k) = 2 * ((-1)^k) / (1 - TGL(k));
end

%last row:
for k = 2:(N-1)
    D1(N,k) = -2 * ((-1)^((N - 1) + k)) / (1 + TGL(k));
end

%first column:
for j = 2:(N-1)
    D1(j,1)= -0.5 * ((-1)^j) / (1-TGL(j));
end

%last column:
for j=2:N-1
    D1(j,N)= 0.5 * ((-1)^((N-1) + j)) / (1 + TGL(j));
end  

%everything else:
for j=2:(N-1)
    for k=2:(N-1)
        if j==k %specific conditions
            D1(j,k) = (-TGL(k)) / (2 * (1 - (TGL(k))^2));
        else
            D1(j,k) = ((-1)^(j + k)) / (TGL(j) - TGL(k));
        end
    end
end
%%%%%END FILLING D1%%%%%

%Second derivative matrix
D2 = D1 * D1;

%our third matrix, which is modified in accordance with BCs
A = eye(N) - (0.5 * dt * alpha * D2);
A_helper = eye(N);
A(1,:) = A_helper(1,:); %modify first and last rows
A(N,:) = A_helper(N,:); %modify first and last rows

%INITIAL CONDITION
U(:,1) = -sin(pi*TGL(:)); %initial condition 
 
 
%Initialises AB for non-linear term
U_dummy1 = U(:,1);
U_dummy2 = D1 * U_dummy1;
Burger1 = U_dummy1.*U_dummy2;
Burger = Burger1;
 
%Initialise Ax=b problem
b = (eye(N) + 0.5 * dt * alpha * D2) * U_dummy1 - dt * 0.5 * (3 * Burger1 - Burger);
b(1) = 0;     % boundary condition 
b(N) = 0;     % boundary condition 
U(:,2) = A\b; %solve Ax=b problem
 
 
%Time loop: includes non-linear term and Ax=b problem
for n=2:N
    U_dummy1 = U(:,n); 
    U_dummy3 = U(:,n-1);
    U_dummy2 = D1*U_dummy1;
    U_dummy4 = D1*U_dummy3;
    Burger1 = U_dummy1.*U_dummy2;
    Burger = U_dummy3.* U_dummy4;
    b=(eye(N) + dt * alpha * D2) * U_dummy1 - dt * 0.5 * (3 * Burger1 - Burger); 
    b(1) = 0;
    b(N) = 0;
    U(:,n+1) = A\b;
end


%PLOTS:
%1: SUPERIMPOSITION
figure
plot(TGL,U(:,1:50)); 
axis([-1 1 -1.5 1.5]);
figure
plot(TGL,U(:,1:58)); 
axis([-1 1 -1.5 1.5]);
figure 
plot(TGL,U(:,58)); 
axis([-1 1 -1.5 1.5]);

% 2: MOVIE
% figure
% for i = 1:2:50
%     plot(TGL,U(:,i)); axis([-1 1 -1.5 1.5]);
%     F(i) = getframe;
% end

% 3: MOVIE TO AVI
% myVideo = VideoWriter('myfile.avi');
% uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
% open(myVideo);
% for i = 1:35
%     plot(TGL,U(:,i)); axis([-1 1 -1.5 1.5])
%     F(i) = getframe;
%     writeVideo(myVideo, F(i));
% end
% close(myVideo);
