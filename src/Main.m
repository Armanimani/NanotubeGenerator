%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmer : Arman Imani 91207404                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is capable of generating nanotubes based on chiral         %
% coordination . The chiral coordinate defines by two vecrots m , n . The %
% input data for m & n should be integer ( so we can roll the graphene    %
% plane ) and also it should be positive. The other two inputs are the    %
% nanotube lentgh and the bond lentgh.                                    %
% The output data includes a graphical output + 2 output files these files%
% are atom_data and bond_data which include the coordinates of atoms and  %
% the bonds between them.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
format long;
%...........................Input Data.....................................
n=input('Enter the Chiral coordinate n>0 & integer : ');
m=input('Enter the Chiral coordinate m>0 & integer : ');
l=input('Enter the nanotube length in nanometer : ');
bl=input('Enter the Bond length in nanometer : ');
a=2*bl*cos(pi/6);
a1=[1;0]*a;
a2=[cos(pi/3);sin(pi/3)]*a;
vector=n*a1+m*a2;
length_vector=sqrt(vector(1,1)^2+vector(2,1)^2);
cos_t=vector(1,1)/length_vector;
sin_t=vector(2,1)/length_vector;
T=[cos_t,sin_t;-sin_t,cos_t];
counter=1;
check=1;
y_temp=-bl/2;
%......................Generating graphene plane...........................
for j=0:5*round(l/bl);
    if (check==3)
        xo=a/2;
        check=check+1;
        y_temp=y_temp+bl/2;
    elseif (check==4)
        x=a/2;
        check=1;
        y_temp=y_temp+bl;
    elseif (check==2)
        xo=0;
        check=check+1;
        y_temp=y_temp+bl;
    else
        xo=0;
        check=check+1;
        y_temp=y_temp+bl/2;
    end
    for i=-5*round(l/a):5*round(l/a);
        x=xo+i*a;
        y=y_temp;
        R=[x;y];
        R_temp=T*R;
        if (and(R_temp(1,1)>=0,R_temp(1,1)<1.01*length_vector))
            if (and(R_temp(2,1)>=0,R_temp(2,1)<l))
                graphene_plane(counter,1)=counter;
                graphene_plane(counter,2)=R_temp(1,1);
                graphene_plane(counter,3)=R_temp(2,1);
                counter=counter+1;
            end
        end
    end
end
%..............Rolling the Graphene plane..........................
Number_atoms=size(graphene_plane,1);
r=length_vector/2/pi;
for i=1:Number_atoms
    teta=graphene_plane(i,2)/r;
    Atom_data(i,1)=i;
    Atom_data(i,2)=r*sin(teta);
    Atom_data(i,3)=graphene_plane(i,3);
    Atom_data(i,4)=r*cos(teta)+r;
end
%.............Generating Bond matrix...............................
for i=1:Number_atoms
    c(1:3,1)=0;    d(1:3,1)=2*a;
    check=1;
    for j=1:Number_atoms
        if (i~=j)
            dx=abs(Atom_data(j,2)-Atom_data(i,2));
            dy=abs(Atom_data(j,3)-Atom_data(i,3));
            dz=abs(Atom_data(j,4)-Atom_data(i,4));
            if (and(and(dx<(1.01*bl),dy<(1.01*bl)),dz<(1.01*bl)))
                dd=sqrt(dx^2+dy^2+dz^2);
                if (dd<max(d))
                    dummy_temp=max(d);
                    dummy=find(d==dummy_temp);
                    dummy=max(dummy);
                    c(dummy,1)=j;
                    d(dummy,1)=dd;
                end
            end
        end
    end
    for k=1:3
        if (c(k)~=0)
            Bond_matrix(i,c(k))=1;
            Bond_matrix(c(k),i)=1;
        end
    end
end
%......................Drawing spheres as atoms..........................
for i=1:Number_atoms
    x=Atom_data(i,2);
    y=Atom_data(i,3);
    z=Atom_data(i,4);
    hold on
    [xx,yy,zz]=sphere();
    rr=0.075*a;
    surf(rr*xx+x,rr*yy+y,rr*zz+z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    colormap([0 0 0.8]);
end
counter=1;
%....................Drawing bonds between atoms...........................
 for i=1:Number_atoms
     x=Atom_data(i,2);  y=Atom_data(i,3);   z=Atom_data(i,4);
     for j=1:i
        if (Bond_matrix(i,j)==1)
            xx=Atom_data(j,2); yy=Atom_data(j,3);  zz=Atom_data(j,4);
            Bond_data(counter,1)=counter;
            Bond_data(counter,2)=i;
            Bond_data(counter,3)=j;
            counter=counter+1;
            line([x,xx],[y,yy],[z,zz],'Color','r','LineWidth',1)
         end
     end
 end
 %.................Generating output files.................................
Atom_data=Atom_data';
Bond_data=Bond_data';
fid1=fopen('Atom_data.txt','wt');
fid2=fopen('Bond_data.txt','wt');
fprintf(fid1,'Atom_number     X(nm)     Y(nm)     Z(nm)     \n');
fprintf(fid2,'Bond_number      Atom1      Atom2     \n');
fprintf(fid1,'%d  %15.10f  %15.10f  %15.10f  \n',Atom_data);
fprintf(fid2,'%d  %d  %d  \n',Bond_data);