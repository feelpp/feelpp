Frequently Asked Questions
==========================

# Execution

## Parallel computing

### How to bind Feel++ mpi processes to cores

For performance reasons, it can be very efficient to bind mpi processes to cores.
You should read the mpirun manual and the binding section. 

Here is the command line to bind the processes:
```
mpirun -np 4 -bind-to-core feelpp_qs_laplacian
```

## Error or Exception messages

### locale::facet::_S_create_c_locale name not valid

On Ubuntu systems and possibly others, locale (language, money, ...) information
might be set and we don't handle that very well at the moment. You might get the
following message at the execution of any Feel++ applications
```
terminate called after throwing an instance of 'std::runtime_error'
  what():  locale::facet::_S_create_c_locale name not valid
```
To fix this, we suggest that you `unset` all environment variables related to
locale information. It includes `LANG`,  `LC_CTYPE` and `LC_ALL`.

Just type
```
unset LANG LC_CTYPE LC_ALL
```
