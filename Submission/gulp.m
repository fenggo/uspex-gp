function status = gulp(generation,Ind_No,Step,commandExecutable)

status = 1;
[a,b] = unix(['mv input_orderd inp_' num2str(generation) '_' num2str(Ind_No) '_' num2str(Step)]);
%[a,b] = unix('mpirun -n 8 gulp<input>output');
[a,b] = unix(commandExecutable);
%[a,b] = unix(['cp gulp.cif optimized_' num2str(generation) '_' num2str(Ind_No) '.cif']);



