%Assignment 1 - COMP558
%QUESTION 2 - EDGE DETECTION
%Author: Camila Alvarez
%Date: 2021 October 3rd


    %Reading Image
    I = imread("skyscrapers.jpg");
    figure('Name','Show Image');
    imshow(I);
    hold off;
    
    I1 = rgb2gray(I); %transform to gray scale img
    figure('Name','RGB to Gray');
    imshow(I1);
    
    [u,v] = size(I1);
    setx = zeros(u,v,3);
    sety = zeros(u,v,3);
    setd = zeros(u,v,3);

    %% 2A) and 2B) Constructing filter
    %Sigma values
    Tx = 500;
    Ty = 10;
    
    %Size of matrix
    A = ((4*Tx) +1);

    %Initializing matrices
    M = zeros(A);
    N = zeros(A);
    O = zeros(A);
    
    %1D Gaussian Formula multiplied by sigma squared to normalize it
    syms x y
    Gaussx =(Tx^2)*(1/(sqrt(2*pi)*Tx))*exp(-(x^2)/(2*(Tx^2)));
    Gaussy =(Ty^2)*(1/(sqrt(2*pi)*Ty))*exp(-(y^2)/(2*(Ty^2))); 
    
    %Computing 2nd partial derivative of Gaussian function w.r.t to x
    G = diff(Gaussx,x,2);
    
    %Initializing rowx
    rowx = {};
    
    %%Populating rowx with 2nd pd of Gaussian function w.r.t x using sigma=Tx
    for i = (-2*Tx):(2*Tx)   
      
            rowx {end+1} = subs(G, {x}, i);    
    end  
    
    %Populating M with Gaussian second derivative values
    for p = 1:A        
            M (:,p)= rowx{p};   
    end

    %Outputting filter for 2A)    
     figure('Name','2A) Filter');
     surfc(M);
     %Outputting 2D height Function for filter
     figure('Name', '2A) 2D height function')
     imagesc(M)
     colormap jet;
     axis image
     axis off;
    
     rowy = {};

    %Populating rowy with 1D Gaussian values using sigma=Ty
      for j = (-2*Tx):(2*Tx)
         rowy {end+1} = (Ty^2)*(1/(sqrt(2*pi)*Ty))*exp(-(j^2)/(2*(Ty^2)));
      end
    
    %Populating N with 1D Gaussian row values
    for p = 1:A        
            N (:,p)= rowy{p};   
    end
    
    Q = N';
    D = Q(:,1);
    
     %Extract columns from M and multiply them by rowy values
        for k = 1:A
            clm = M(:,k);   
            for r = 1:A
                temp = clm(r)*D(r);
                clm(r) = temp;
            end
            O(:,k)=clm;
        end
      %Outputting 2B) Filter    
      figure('Name','2B) Filter');
      surfc(O);

    
        %% 2C) ROTATION
        B = ((8*Tx) +1);
    
        % Creating bigger filter for rotation purposes
        rotation = zeros(B); %bigger size filter that will be rotated
        E = zeros(B); %1D Gaussian tapering filter
    
        %Initializing rowz and rowg
        rowz = {};
        rowg = {};
    
       %Populating rowx with 2nd pd of Gaussian function w.r.t x using sigma=Tx
        for i = (-4*Tx):(4*Tx)   
            rowz {end+1} = subs(G, {x}, i);    
        end  
    
        % Populating rotation with Gaussian second derivative values
        for p = 1:B        
            rotation (:,p) = rowz {p};   
        end
    
        % Populating rowy with 1D Gaussian values using sigma=Ty
        for j = (-4*Tx):(4*Tx)
            rowg {end+1} = (Ty^2)*(1/(sqrt(2*pi)*Ty))*exp(-(j^2)/(2*(Ty^2)));
        end
    
        % Populating E with 1D Gaussian row values
        for p = 1:B        
            E(:,p)= rowg{p};   
        end
    
        Tr = E';
        F = Tr(:,1);
    
        % Extract columns from rotation and taper them 
        for k = 1:B
            clm = rotation(:,k);  
            for r = 1:B
                temp = clm(r)*F(r);
                clm(r) = temp;
            end
            rotation(:,k)=clm;
        end
    
        % ROTATION PROCEDURE to generate rotated filters by 0, 45, 90 and 135 degrees 
        % Can be modified to work with any angles
          theta_array = [0 (pi/2) (pi/4) (3*(pi/4))]; 
          centerS = (A+1)/2;
        
          for n = 1:length(theta_array) 
               countery = (4*Tx) +1;
               theta = theta_array(n);
               rotated = zeros(A); %same size as original matrix from 2
                  for i= 1:B %Iterating through rows (y axis)        
                    countery = countery -  1;
                    counterx = (-4*Tx) - 1; 
                        for j = 1:B %Iterating thorugh columns (x axis)        
                            counterx = counterx + 1;
                            rx = round(counterx*cos(theta) - countery*sin(theta));
                            ry = round(counterx*sin(theta) + countery*cos(theta));                
                            ax =centerS+rx;
                            ay = centerS+ry;                
                            if (ry >= (-2*Tx) ) && (ry <=(2*Tx))
                                 if (rx >= (-2*Tx)) && (rx <= (2*Tx))      
                                    rotated(ay,ax) = rotation (i,j);
                                 end             
                            end
                        end
                  end
            
            
              %2D) OUTPUT 
              %Filtering image
              Iconvx = conv2(I1,rotated, 'same');
              
              %Outputting results of filter and filtered image
              strg = append('2) Rotated Filter and Filtered Image for ', num2str(theta),' angle in radians.');
              figure('Name',strg);
              subplot(1,2,1),surfc(rotated);
              subplot(1,2,2),imshow(Iconvx);  
              

              %cirshift in x-direction
              %shift whole matrix to the left by 1
              setx (:,:,1) = (sign(circshift(Iconvx,-1,1))~=sign(Iconvx))*0.5;
    
              %cirshift with y-direction
              %shift whole matrix up by 1
              sety (:,:,2) = (sign(circshift(Iconvx,-1,2))~=sign(Iconvx))*0.5;
    
              %cirshift with diagonal-direction
              %shift whole matrix up by 1
              setd (:,:,3) = (sign(circshift(Iconvx,[-1,-1]))~=sign(Iconvx))*0.5; 

              P = imfuse(setx (:,:,1),sety (:,:,2));
              S = imfuse(P,setd(:,:,3));

              strg2 = append('2C)Fused zero crossings for ', num2str(theta),' angle in radians.');
              figure('Name',strg2);
              imshow(S)


          end
  
            
    
    
          %% 2E) LoG    
          LoG = fspecial('log',40,6);
          %Since LoG is rotationally invariant LoG will be the same whenever it´s rotated
           
          Icovx = conv2(I1,LoG, 'same');
          %Initializing binary image
              cross = zeros(size(I1));
              %cirshift in x-direction   
              %shift whole matrix to the left by 1
              cross(:,:,1) = (sign(circshift(Icovx,-1,1))~=sign(Icovx));
 
              %cirshift with y-direction
              %shift whole matrix up by 1
              cross(:,:,2) = (sign(circshift(Icovx,-1,2))~=sign(Icovx));
    
              %cirshift with diagonal-direction
              %shift whole matrix up by 1
              cross(:,:,3) = sign(circshift(Icovx,[-1,-1]))~=sign(Icovx);

              figure('Name','2D) Zero-Crossings for Laplace of Gaussian');
              subplot(2,2,1),imshow(cross(:,:,1));
              subplot(2,2,2),imshow(cross(:,:,2));
              subplot(2,2,3),imshow(cross(:,:,3));       

