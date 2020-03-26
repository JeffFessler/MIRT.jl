A student used the following steps
on a Windows 64-bit machine
to set up multi-thread processing.

1. Click on Windows key on your keyboard.
Type "environment variables" that will take you to
"system environment variables" dialogue box.

2. Add new, `variable name= JULIA_NUM_THREADS, variable value= %NUMBER_OF_PROCESSORS%` 
This automatically gets converted to 8 after I hit Enter.

3. Open Windows Powershell (cmd),
go to your Julia bin directory using `cd` command.
For me, this was at
`C:\Users\Username\AppData\Local\Julia-1.3.1\bin`

4. Type `.\julia` (or `julia` - whichever your cmd trusts to open).

5. A Julia terminal will open within Windows cmd.
Typing `Threads.nthreads()`
should give you the same value as the environment variable you set above.

Further experiments: 

1. Through Julia command-line
`Threads.nthreads()` returns 8,
while in Atom,
`Threads.nthreads()` returns 4. 
This seems to be the effect of using virtual versus physically available cores.

2. Two more commands might fetch your attention:
`using Distributednworkers()`
This value is by default 1,
and can be changed in Windows terminal while starting Julia:
`.\julia --banner=no -p 3`

This deploys 3 workers in Julia. 
