Nkx = length(kPoints(1,:));
Nky = length(kPoints(2,:));

j0 = Nkx/2 +1 ;
i0 = Nky/2 + 1;

thresh = .025;

sig = .05;
mu = 0;

FEig = exp(-((MyEig-mu).^2)/(sig^2)); %check if exp is entry

track = 0;
for i = 1:50
        ind = find(FEig(i,:)>thresh);
        track = max(track, length(ind));
end

[B, I] = sort(FEig(:, j0));
ind = I;%((end-track +1):end);
FE = diag(FEig(ind, j0)); %Track largest values of FEig instead of all, right size
track = length(ind);


Psi = MyEigVecs(:, ind, j0);%, i0);

[Q, R, p] = sRRQR_rank(FE*Psi', 1, 3);

% Constructs nb x nb matrix Psi

% Psitild = zeros(M*N, track,length(ky),length(kx));
% Psi = zeros(M*N, track,length(ky),length(kx));

Psitild = zeros(13468,50,length(Nkx));
%relabel indices so j i consistent with elsewhere
%Single loop
    for j = 1:length(Nkx)
        [~, I] = sort(FEig(:, j));
        ind = I((end-track +1):end);
        Psi = MyEigVecs(:, ind, j );
        FE = diag(FEig(ind, j));
        E = FE*Psi(p(1:track),:)';
        B = sqrtm(E'*E);
        U = E/B;
        Psitild(:,:,j) = Psi*U;
    end


Wann = ifft(ifft(Psitild, [], 3), [], 4);
Wann = fftshift(Wann, 3);
Wann = fftshift(Wann, 4);

W = zeros(length(Nky)^2,length(Nkx)^2, track);

for n = 1:Nky
   for m = 1:Nkx
       for b = 1:track
          for i = 1:length(Nky)
              for j = 1:length(Nkx)
              W(n + length(Nky)*(j-1), m + length(Nkx)*(i-1), b) = Wann(n + length(Nky)*(m-1), b,j,i);
              end
          end    
       end
   end
end

for i = 1:1%track
    figure(i)
    clf()
    surf(abs(W(:,:,i)))
end