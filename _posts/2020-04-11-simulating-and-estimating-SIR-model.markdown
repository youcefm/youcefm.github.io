---
layout: post
title:  "Bayesian Simulation of a SIR Model in Python"
date:   2020-04-11 
categories: bayesian statistics| dynamic models| epidimiology
---


### Dynamic System Equations
Lets start by writing down the functions that describe how the population evolves as the virus spreads.

{% highlight python %}

def dy_dt(y, z, beta, sigma):
    "change in proportion of population infectious"
    y = min(max(0,y),1)
    z = min(max(0,z),1)
    beta = max(0,beta)
    sigma = max(0,sigma)
    
    return beta*y*(1-z) - sigma*y
{% endhighlight %}

Check out the [Jekyll docs][jekyll-docs] for more info on how to get the most out of Jekyll. File all bugs/feature requests at [Jekyll’s GitHub repo][jekyll-gh]. If you have questions, you can ask them on [Jekyll Talk][jekyll-talk].

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/
