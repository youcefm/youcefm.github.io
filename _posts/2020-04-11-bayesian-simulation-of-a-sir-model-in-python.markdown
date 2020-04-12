---
layout: post
title:  "Bayesian Simulation of a SIR Model in Python"
date:   2020-04-11 
categories: dynamic-models-1 epidimiology-2 differential-equations-3
mathjax: true
---

In this post, I want to explore the canonical epidimioligical model of virus spread known as the Susceptible-Infectious-Recovered (SIR) model. I draw on papers that have been published in the last month on COVID-19 to build a full Bayesian simulation of a SIR model fitted to the daily data on deaths, using Python. 

Topics: Epidimiology, Dynamic Models, Bayesian Statistics, MCMC, Python.

## Introduction
In March 2020, we saw a number of papers published on [Medrxiv](https://www.medrxiv.org/) trying to understand the spread of the virus SARS-COV2, the coronavirus responsible of COVID-19 respiratory illness. Many of these papers use some type ofa Susceptible-Infectious-Recovered framework to model the disease and analyze data.

In a series of 3 posts I will walk you through the process of defining a SIR model to simulate the spread of a virus in a population, estimate the parameters of the model using Bayesian methods, and finally making predictions from the model using bayesian predictive posteriors. All of this is done in Python. Existing academic models of this kind are usually built in C. This is because the MCMC simulation loop can be implemented much more efficiently in C than Python. But, building such models in Python can allow much faster iteration speed, and allow the exploration of amny more model variations.

This post, the first in the series, walks you through the basic ingredients of SIR models and shows you how to build a variant of the model defined in [this paper](https://www.medrxiv.org/content/10.1101/2020.03.24.20042291v1.full.pdf).

In the second post, I introduce data. I explain how to define priors on the model parameters as well as a liklihood function for the data given the parameters. I then build an Markov Chain Monte Carlo (MCMC) algorithm from scratch to estimate the posetrior distribution of the parameters given the model, the priors and the data.

In the last post, I show you how to use the output of the posterior draws from the MCMC simulation to then predict how the virus spread will evolve in the future, building a full range of scenarios taking parameter uncertainty into account.

### The Basic Ingredients of SIR Models

The first proposed SIR model seems to date back to a [1927 paper by W. O. Kermack and A. G. McKendrick](https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1927.0118). For a simple step by step explanation of SIR models, I encourage you to [read this series of posts](https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model) from the Mathematical Association of America. Here, I will explain the minimum needed to implement our model in Python.

A SIR model is the simplest way to mathematically describe the spread of a virus in a population over time. The model treats the population as a homogenous blob with random interactions between individuals, abstracting away the specific topology of the network underlying social interactions.

It starts by considering a population susceptible to be infected by the virus (say everyone in a give city or country). Once the virus is introduced in the population (first infection), infected individuals transmit the virus to other susceptible individuals at the transmission rate $\beta$ per unit of time. Infected Individuals then "recover" from the infection at a rate $\sigma$ per unit of time, implying an infectious period of $\div{1}/{\sigma}$. This can be jarring to realize, but the term recovery in these models simply means the individual is no longer infectious; they could still feel terrible and eventually die.

A critical quantity in SIR models is the base reproduction rate of the virus, $R_0$. It is the number of new infections resulting from every existing infection. This quantity is equal to the transmission rate times the infectious period, $R_0 = \div{\beta}/{\sigma}$ . Intuitively, the virus reproduces more if it is more likely to transmit from an infected individual to susceptible individuals (higher $\beta$), and it also reproduces more if infected individuals remain in the infectious pool for longer periods of time. 

The number of deaths from the virus can be obtained in the SIR model by adding an extra parameter: the Infection Fatality Rate (IFR). This is the number of infected individuals that "recover" by dying from the illness the virus causes.

Now that we understand what SIR models are made of, lets build one using a system of differential equations. The goal is to define equations that, given initial conditions, can simulate the full evoluation of the virus spread in the population.

### Differential Equations Describing a SIR Model
Differential equations are very useful for describing how some quantity changes over time. They can be used to describe how the virus spreads over time in the context of a SIR model. As discussed in the previous section, the SIR model is composed of three pools: the fraction of the population in the infectious pool \\(y(t)\\), the fraction of the population in the recovered pool \\(z(t)\\) (no longer infectious), and the fraction of the population susceptible ( \\(s(t) = 1 -z(t) -y(t)\\) ). 
To fully describe this system, we need equations that tell us how the infectious and recovered populations change at any point in time.  

$$\div{dy}{dt} = \beta*y*s - \sigma*y$$

$$\div{dz}{dt} = \sigma*y$$

$$\div{ds}{dt} = -\beta*y*s $$

Now that we have these equations, we can define an equation for the number of deaths from the virus. We need to define parameters that determine how infected individuals can die. $/rho$ is the rate of severe disease from infection, and $\theta$ is the risk of death conditional on severe disease. The product of these two parameters is the IFR. Given a total population N, the equation for cumulative deaths over time is:

$\Omega(t) = N*\rho*\theta*z_{t - \psy}$

The initial conditions are given by $s(0) = 1$ , $y(0) = 1/N$ , $z(0) = 0$.

These equations are simulated in Python with the following code

{% highlight python %}

def dy_dt(y, s, beta, sigma):
    "change in proportion of population infectious"
    y, s, beta, sigma = min(max(0,y),1), min(max(0,s),1), max(0,beta), max(0,sigma)
    return beta*y*s - sigma*y
    
def dz_dt(y, sigma):
    "change in proportion of population no longer infectious"
    y, sigma = min(max(0,y),1), max(0,sigma)
    return sigma*y
    
def ds_dt(y, s, beta):
    "change in proportion of population susceptible"
    y, s, beta = min(max(0,y),1), min(max(0,s),1), max(0,beta)
    return -beta*y*s

def cum_deaths(N, rho, theta, z_lag):
    "cumulative deaths"
    rho, theta, z_lag = min(max(0,rho),1), min(max(0,theta),1), min(max(0,z_lag),1)
    return np.round(N*rho*theta*z_lag)

def simulate_data(inverse_sigma=4.5, theta=0.14, rho=0.01, R=2.75, psy=17):
    "simulate output of SIR model"
    
    #constants
    N=66*10**6          ## population size
    time_span=90        ## simulation period
    first_case = 10     ## time of introduction of first infection
        
    #parameters
    sigma = 1/inverse_sigma
    beta = sigma*R 
    lag = max(0, int(np.round(psy)))
    y, z, s, deaths = [0], [0], [1], [0]
    
    for i in range(0, time_span):
        if i <first_case:
            y.append(0)
            z.append(0)
            s.append(1)
            deaths.append(0)
        elif i == first_case:
            y.append(1./N)
            z.append(0)
            s.append(1)
            deaths.append(0)
        else:
            y.append(y[i] + dy_dt(y[i], s[i], beta, sigma) )
            z.append(z[i] + dz_dt(y[i], sigma) )
            s.append(s[i] + ds_dt(y[i], s[i], beta) )
            deaths.append(cum_deaths(N, rho, theta, z[max(0, i - lag)]) ) ## works because z is always 0 at index 0

        return [deaths, y, z, s]
{% endhighlight %}

### Solution Curves

{% highlight python %}

model_output = simulate_data()
df = pd.DataFrame(zip(model_output[0], model_output[1], model_output[2], model_output[3], range(0,100)), 
columns=['cum_deaths', 'prop_infectious', 'prop_recovered', 'prop_susceptible', 'time'])
df.plot(x='time', y=['prop_infectious', 'prop_recovered', 'prop_susceptible'], 
title= 'Solution Curves of SIR Model')

{% endhighlight %}

[GRAPH HERE]

### Compare Model To Real Data from the UK

Now lets compare our model output for number of death to real data from the UK. I downloaded the data from this [John Hopkins University Github repo](https://github.com/CSSEGISandData/COVID-19). Here I use the archived data as of March 25 2020. Feel free to replace by more up to date data.

{% highlight python %}

## prep UK data
data = pd.read_csv('covid_deaths.csv')
uk_data = data[(data['Country/Region'] == 'United Kingdom') & (data['Province/State'] == 'United Kingdom')]

for column in ['Province/State', 'Country/Region', 'Lat', 'Long']:
    del uk_data[column]

uk_data = uk_data.melt(var_name ='date', value_name = 'deaths')
uk_data['cum_deaths'] = uk_data.deaths.cumsum()
uk_data = uk_data = uk_data.loc[0:57]

{% endhighlight %}

Now lets compare the actual UK deaths to the model output. The model output will depend on the specific parameter values chosen. I calibrated the parameters to relflect the UK situation reasonably well.

{% highlight python %}

df_compare = pd.DataFrame(zip(range(0, 97), simulate_data()[0], uk_data.cum_deaths), columns=['time', 'model_output', 'data'])

import matplotlib.pyplot as plt
model, = plt.plot(df_compare.time, df_compare.model_output, 'k^', label='model output')
data, = plt.plot( df_compare.time, df_compare.data, 'r+', label='observed data')
plt.legend(handles=[model, data])
plt.title('Actual Deaths vs Model Output')
plt.xlabel('time')
plt.ylabel('cumulative deaths')

{% endhighlight %}

[GRAPH HERE]
