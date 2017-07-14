
%Please Pick your favourite ALPHA (i.e., diffusion coefficient) and RUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.0001;%diffusion coefficient%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


J = 100;    %space steps
N = 100;    %time steps

dx = 2 / (J-1); %spatial resolution
x = -1:dx:1;    %x range

%CFL conditions:
if alpha <0.0001
    t_max = 1;
elseif alpha < 0.001
    t_max = 1.2;
elseif alpha < 0.01
    t_max = 3;
elseif alpha < 0.1
    t_max = 5;
elseif alpha < 2
    t_max = 1;    
end
 
dt = t_max / N; %space step size

  

U = zeros(N,J);         %Initialise matrix for data storage
U(1,:) = -sin(pi*x(:)); %Initial condition is Sin(pi*x) s.t. x \in [-1 1]
  
%Initalise function using lower order method
n=1;
%go to fourier space with IC
FF = fft(U(1,:));                   
%likewise for non-linear term and get next step (THIS IS ADAMS-BASHFORTH)
dummy_1 = FF(n,1:J);
dummy_2 = 1i*[0:J/2-1 0 -J/2+1:-1].*dummy_1;
dummy_1 = real(ifft(dummy_1));
dummy_2 = real(ifft(dummy_2));
Burgers_current = dummy_1.*dummy_2;
Burgers_current = fft(Burgers_current);

%THIS TESTS HEAT EQUATION ONLY:
%Burgers_current = zeros(1,J);   

for k=1:J
    U(n+1,k)=((1+(dt*alpha*(k^2)))^(-1))*( FF(n,k) - dt * Burgers_Current(k));
end

%Shift values for next timestep
Burgers_previous=Burgers_current;
  
%Begin Time loop
for n=2:N
    FF=U(n,:);
    
    dummy_3 = U(n,1:J);
    dummy_4 = 1i*[0:J/2-1 0 -J/2+1:-1].*dummy_3;
    dummy_3 = real(ifft(dummy_3));
    dummy_4 = real(ifft(dummy_4));
    Burger_current = dummy_3.*dummy_4;
    Burgers_current = fft(Burger_current);
 
    %THIS TESTS HEAT EQUATION ONLY:
    %Burgers_current = zeros(1,J);
    
    for k= 1:J    
      U(n+1,k) = (1/(1-(dt*alpha*(-k^2))*.5))*(FF(k)*(1+(dt*alpha*(-k^2))*.5)...
        -(dt*.5)*(3*Burgers_current(k)-Burgers_previous(k))); 
    end
    Burgers_previous = Burgers_current; %shift values for next time step
end
  
%pull that which is in Fourier space out of it (first step is alread done)
for j=2:N
    U(j,:)=real(ifft(U(j,:)));
end
  
U=real(U); %just for insurance...matlab is often a failure...should really be in Julia






%plots:
%1: SUPERIMPOSITION
% for i = 1:1:N
%     plot(x,U(i,:)); hold on;
% end

%2: MOVIE
figure
for i = 1:N
    plot(x,U(i,:)); axis([-1 1 -1 1]);
    F(i) = getframe;
end

%3: MOVIE TO AVI
% myVideo = VideoWriter('myfile.avi');
% uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
% open(myVideo);
% for i = 1:N
%     plot(x,U(i,:)); axis([-1 1 -1.5 1.5])
%     F(i) = getframe;
%     writeVideo(myVideo, F(i));
% end
% close(myVideo);
  