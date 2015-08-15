%Promena parametra frekvencije

c = 3*10^8; %brzina svetlosti [m/s]
lambda=998.85e-9;
w=2*pi*c/lambda; %ugaona frekvencija
wmin=0.8*w;
wmax=1.2*w;
v = [wmin:0.001*w:wmax]; %posmatrani ugaono-frekvencijski domen

a = 156e-9; %duzina perioda
Nperiod = 20; %broj perioda u 1D fotonskom kristalu do defekta i od defekta
n1 = 2.9; %indeks prelamanja prvog materijala
d1 = 70e-9; %duzina prvog materijala u periodu
n2 = 3.57; %indeks prelamanja drugog materijala
d2 = a-d1; %duzina drugog materijala u periodu
n_0 = 1; %indeks prelamanja okoline
%  d=4*d1; %defekt u strukturi
d=2*d1;
% dmin=d1;
% dmax=4*d1;
% d=[dmin:0.01*d1:dmax];
step=0.01;

step=step*a;
I=sqrt(-1);




for k=1:length(v);
    
C_1=[d1 d2]; %do defekta 
C_2=[d2 d1]; %od defekta 
D_1=repmat(C_1,[1 Nperiod]); %niz sa vrednostima duzina do defekta
D_2=repmat(C_2,[1 Nperiod]); %niz sa vrednostima duzina od defekta
L=[D_1 d D_2]; %niz sa vrednostima duzina u prostoru sa defektom

rmin=0;
rmax=sum(L);
r=[rmin:step:rmax];


n_m=[n1 n2];
n=repmat(n_m,[1 Nperiod]);
n=[n n1 fliplr(n)]; %niz sa vrednostima indeksa prelamanja

%Folder


%"interface" matrix, matrica prelaza izmedju materijala 1 i okoline 
A = 1/2*[1+n_0/n1 1-n_0/n1;1-n_0/n1 1+n_0/n1];
A_prim=1/2*[1+n1/n_0 1-n1/n_0;1-n1/n_0 1+n1/n_0];
 
%"interface" matrix, matrica prelaza izmedju 2 sredine
W_1=1/2*[1+n1/n2 1-n1/n2;1-n1/n2 1+n1/n2];
W_2=1/2*[1+n2/n1 1-n2/n1;1-n2/n1 1+n2/n1];



broj_slojeva=size(L,2); %ukupan broj slojeva u kristalu
L1=L;
broj_tacaka=length(r); %broj tacaka u kojima zelim da nadjem vrednost polja
 
S=zeros(2,broj_tacaka); %matrica stanja elektricnog i magnetnog polja u krajnjoj tacki
 

 S(:,broj_tacaka)=inv([1 1; -n1 n1])*[1 1; -n_0 n_0]*[1;0]; %stanje u tacki na kraju strukture    
 j=broj_tacaka-1; %za prvi step

 B=eye(2);
 cuc=1;
 
 kumulativnoL=cumsum(L);
 	
 curr_cell=broj_slojeva;
 
 mkdir(strcat('Plotovi'));
 

 
   for i=(broj_tacaka-1):(-1):1
    
    if(mod(i,1000))
        i/broj_tacaka
    end
   
    
    j=find(kumulativnoL>=r(i), 1 );
    
    if( j< curr_cell)
        
        if (mod(curr_cell,2)==1) 
            
        B=W_2*get_transfer_matrix1(n(curr_cell),v(k),-L(curr_cell),c)*B;
        
        else
            B=W_1*get_transfer_matrix1(n(curr_cell),v(k),-L(curr_cell),c)*B;
        end
        curr_cell=curr_cell-1;
        
    end
    
    S(:,i)=get_transfer_matrix1(n(curr_cell),v(k), r(i) - kumulativnoL(j) ,c)*B*S(:,broj_tacaka);
   

end

 


h=figure(k);
grid on;
Q=(abs(S(1,:)+S(2,:)).^2);
najm=min(Q);
najv=max(Q);
plot(r,Q,'r');
% plot(r(14001:27795),Q(14001:27795),'r');
 %plot(r(14001:26898),Q(14001:26898),'r');

xlabel('x[m]')
ylabel('|E|^2');
title(['Jacina elektricnog polja ' , ', d = ' num2str(d),', w = ' num2str(v(k))]);

% for i=1:broj_slojeva
%     line([kumulativnoL(i), kumulativnoL(i)], [najm najv])
% end


 saveas(h,sprintf('figure_%d.jpg',k));
 close(h)

 end
 
 