---
layout: post
title:  "Bayesian Simulation of a SIR Model in Python"
date:   2020-04-11 
categories: bayesian-statistics dynamic-models epidimiology
---

In this post, I want to explore the canonical epidimioligical model of virus spread known as the Susceptible-Infectious-Recovered (SIR) model. I draw on papers that have been published in the last month on COVID-19 to build a full Bayesian simulation of a SIR model fitted to the daily data on deaths, using Python.

## Introduction
In March 2020, we saw a number of papers published on Medrxiv](https://www.medrxiv.org/) trying to understand the spread of the virus SARS-COV2, the coronavirus responsible of COVID-19 respiratory illness.

### Dynamic System of Equations Describing a SIR Model
Lets start by writing down the functions that describe how the population evolves as the virus spreads.
$ \div{dy}{dt} = \beta*(1-z) - \sigma*y $

{% highlight python %}

def dy_dt(y, z, beta, sigma):
    "change in proportion of population infectious"
    y = min(max(0,y),1)
    z = min(max(0,z),1)
    beta = max(0,beta)
    sigma = max(0,sigma)
    
    return beta*y*(1-z) - sigma*y
def dz_dt(y, z, beta):
    "change in proportion of population no longer infectious"
    y = min(max(0,y),1)
    z = min(max(0,z),1)
    beta = max(0,beta)
    
    return beta*y*(1-z)

def cum_deaths(N, rho, theta, z_lag):
    "cumulative deaths"
    rho = min(max(0,rho),1)
    theta = min(max(0,theta),1)
    z_lag = min(max(0,z_lag),1)
    
    return np.round(N*rho*theta*z_lag)

def simulate_data(inverse_sigma=4.5, theta=0.14, rho=0.01, R=2.75, psy=17):
"Simulate Output of SIR Model"
        # constants
        N=66*10**6
        time_span=57         ## needs more thought
        first_case = 2       ## needs more thought
        #import pdb; pdb.set_trace()
        
        #parameters
        sigma = 1/inverse_sigma
        beta = sigma*R 
        lag = max(0, int(np.round(psy)))
        y=[0]
        z = [0]
        deaths = [0]
        for i in range(0, time_span):
            if i <first_case -1:
                y.append(0)
                z.append(0)
                deaths.append(0)
            elif i == first_case -1:
                y.append(1./(N))
                z.append(0)
                deaths.append(0)
            else:
                y.append(y[i] + dy_dt(y[i], z[i], beta, sigma))
                z.append(z[i] + dz_dt(y[i], z[i], beta))
                deaths.append(cum_deaths(N, rho, theta, z[max(0, i - lag)])) ## works because z is always 0 at index 0
        
        return [deaths, y, z]
{% endhighlight %}

Check out the [Jekyll docs][jekyll-docs] for more info on how to get the most out of Jekyll. File all bugs/feature requests at [Jekyll’s GitHub repo][jekyll-gh]. If you have questions, you can ask them on [Jekyll Talk][jekyll-talk].

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/
