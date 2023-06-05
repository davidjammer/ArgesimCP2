function [u,t] = run(location, type, n, Nx, Ny, REP, saveFlag)
%%Description
%run lbm local or on Seneca in seq. or parallel
%%
%location: "local" or "cluster"
%type: "seq" or "para"
%n: number of worker
%Nx and Ny: size
%REP: number of repetitions
%saveFlag: true -> u is stored in u.dat
%%
	if(location == "local")
		
        if(type == "seq")
            [u,t]=lbm_seq(Nx,Ny,REP);
        else
            fprintf("run local with %d Worker(s)\n", n); 
            pool = parpool("local",n);
            [u,t]=lbm_para1(Nx,Ny,REP);
            delete(pool);
        end
	elseif(location == "cluster")
	    fprintf("run on cluster with %d Worker(s)\n", n);
	    cluster = parcluster('Seneca');
        if(type == "seq")
            job = batch(cluster,	@lbm_seq, 2, {Nx, Ny, REP},'Pool', 1, ...
						'AutoAddClientPath',false);
        else
            job = batch(cluster,	@lbm_para1, 2, {Nx, Ny, REP},'Pool', n, ...
						'AttachedFiles', 'shift_N' ,...
						'AttachedFiles', 'shift_E' ,...
						'AttachedFiles', 'shift_S' ,...
						'AttachedFiles', 'shift_W' ,...
                        'AttachedFiles', 'shift_NE' ,...
						'AttachedFiles', 'shift_SE' ,...
						'AttachedFiles', 'shift_SW' ,...
						'AttachedFiles', 'shift_NW' ,...
                        'AttachedFiles', 'cart2D_shift' ,...
						'AutoAddClientPath',false);
        end
	    
		wait(job);
		x = fetchOutputs(job);
		u = x{1};
        t = x{2};
		
		delete(job);
	end
	fprintf("Elapsed time is %f seconds\n", t);
	
	plot_u(u, REP); 
    if nargin == 7 && saveFlag
      save("u.dat", "u", "-ascii");
    end
end

