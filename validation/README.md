# Validation notes

## Benchmarking

Expected values for the supplied synthetic case:

- initial crossflow: 799.9205 sm3/d
- equilibrium R1 pressure: 256.9231 bar
- equilibrium R2 pressure: 236.9231 bar
- equilibrium transferred volume: 8,699.80 m3
- time constant: 10.8758 d
- time to 0.01 bar excess potential: 75.1275 d

## Eclipse/OPM Flow comparison

Observed in the validated synthetic deck:

- initial connection water flow: about 799.50 sm3/d
- late-time crossflow approaches zero
- final R1 pressure: about 256.92 bar
- final R2 pressure: about 236.92 bar
- regional WIP redistribution: about 8,703–8,704 sm3
- inter-region grid cumulative flow: zero

### Timestep sensitivity

A coarse schedule with hourly report steps for the first day and then daily steps produced approximately:

- well cumulative outflow: 8,010 sm3
- material-balance error: 693–694 sm3

Using:

```text
TSTEP
  1000*0.1
/
```

reduced the reported values to approximately:

- well cumulative outflow: 8,624 sm3
- material-balance error: 80 sm3

The total remains approximately:

```text
8624 + 80 = 8704 sm3
```



