%Raspodela E polja 4; 

%Raspodela elektricnog polja, granicne povrsine
a = 1; %duzina perioda
Nperiod = 10; %broj perioda u Bragovom reflektoru
n1 = 2.5; %indeks prelamanja prvog materijala
d1 = 0.8* a; %duzina prvog materijala u periodu
n2 = 1; %indeks prelamanja drugog materijala
d2 = (1-d1)*a; %duzina drugog materijala u periodu
n = 1; %indeks prelamanja okoline

 c = 3*10^8; %brzina svetlosti [m/s]
%  wmin = 0;
%  wmax = 0.4*2*pi*c;
%  w = [wmin:10^6:wmax]; %posmatrani ugaono-frekvencijski domen

I = 1i; %imaginaran broj
t = [];

%Elektromagnetni talas
w=3*2*pi*c; %ugaona frekvencija
rmin=0.1;
rmax=Nperiod*a;
r=[rmin:0.1:rmax];


   %matrica za jedan period strukture
   Tperiod=[cos(n2*w*d2/c), I*sin(n2*w*d2/c)/ n2; I*sin(n2*w*d2/c)*n2, cos(n2*w*d2/c)]*[cos(n1*w*d1/c), I*sin(n1*w*d1/c)/ n1; I*sin(n1*w*d1/c)*n1, cos(n1*w*d1/c)];   
   Tn = Tperiod;
    for j = 1:Nperiod-1
        Tn = Tn*Tperiod; %matrica za celu strukturu
    end
    
    T_x = Tn; %za vrednost x=0;
    
     A = [1, 1;- n, n];
           M=T_x*A;   
           Minv = [1 0; 0 1]/M;
           P=[1 0]*Minv;
           F=P(1); %napred propagirajuci talas
           B=P(2); %nazad propagirajuci talas
           
            l=1; %brojac elemenata u nizu E koji sadrzi vrednosti jacina elektricnog polja u svakoj tacki x
            E_x=B; %za x=0
            E(l)=E_x;
    
    l=2;
    j=0; %brojac perioda 
    
 for x=r
     
        if x>0 && x<(j+1)*d1+j*d2;
           T=Tperiod^(Nperiod-j); %matrica za elektromagnetni talas u Nperiod-j perioda do granicne povrsine
           T_2=[cos(n2*w*d2/c), I*sin(n2*w*d2/c)/ n1; I*sin(n2*w*d2/c)*n1, cos(n2*w*d2/c)]; %matrica za sloj materijala indeksa prelamanja n2
           T_y1=[cos(n1*w*(d1-x)/c), I*sin(n1*w*(d1-x)/c)/ n1; I*sin(n1*w*(d1-x)/c)*n1, cos(n1*w*(d1-x)/c)];  %koordinata x u prvom materijalu, prostiranje talasa
           T_x=T_y1*T_2*T; %matrica od tacke x do kraja strukture
           A = [1, 1;- n, n];
           M=T_x*A;   
           Minv = [1 0; 0 1]/M;
           P=[1 0]*Minv;
           F=P(1); %napred propagirajuci talas
           B=P(2); %nazad propagirajuci talas
            
            E_x=F*(cos(n1*w*x/c)-I*sin(n1*w*x/c))+B*(cos(n1*w*x/c)+I*sin(n1*w*x/c));
   
        else  T_y2=[cos(n2*w*(a-x)/c), I*sin(n2*w*(a-x)/c)/ n2; I*sin(n2*w*(a-x)/c)*n2, cos(n2*w*(a-x)/c)];
            T=Tperiod^(Nperiod-j);
            T_x=T_y2*T; %matrica od tacke x do kraja strukture,ukoliko je koordinata x u drugom materijalu
            A = [1, 1;- n, n];
            M=T_x*A;   
            Minv = [1 0; 0 1]/M;
            P=[1 0]*Minv;
            F=P(1);
            B=P(2); 
            
            E_x=F*(cos(n2*w*x/c)-I*sin(n2*w*x/c))+B*(cos(n2*w*x/c)+I*sin(n2*w*x/c));
        end
        l=2;
        E(2)=E_x;
             
          if  x<(j+1)*a 
             j=0;
          else j=j+1;
          end
l=l+1;
 
 end
 
 
figure()
hold on;
grid on;
plot(r,abs(E).^2,'r','linewidth',2);
xlabel('Rastojanje [rad/s]')
ylabel('Jacina elektricnog polja');
title(['Jacina elektricnog polja Np = ' int2str(Nperiod) ', n1 = ' num2str(n1) ', n2 = ' num2str(n2) ', n = ' num2str(n), ', d1 = ' num2str(d1)])
hold off;