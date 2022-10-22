using Plots, FFTW, DSP

# Q1
function procedure_1()
    function Xt(t, W)
        return sin(W * t) / pi * t
    end

    function Ut(t)
        if t >0.0
            return 1.0
        else
            return 0.0
        end
    end
    n = 2001
    fs = 1/100
    wc = 100.0
    a = 0.75

    w = -10.0:fs:10.0
    t = -10.0:fs:10.0

    xt = Xt(t,3)


    ut = U.(t)

    function Q1(xt)
        ht = extp.(-a * t) .* ut

        y1 = xt .* cos.(wc * t)
        y2 = xt .* sin.(wc * t)
        w1 = conv(y1, ht)[:2001]
        w2 = conv(y2, ht)[:2001]

        xt1 = w1 .* cos.(wc * t)
        xt2 = w2 .* sin.(wc * t)
        return y1, y2, w1, w2, xt1 + xt2
    end

    y1t, y2t, w1t, w2t, yt = Q1(xt)

    Y1w = fft(y1t)
    Y2w = fft(y2t)
    W1w = fft(w1t)
    W2w = fft(w2t)
    Yw = fft(yt)

    # Plotting
    #plot!(t, real(q1t), label="q(t)", ylim = [-0.0001, 0.0001], xlim = [-1, 1])
end

## Q2

function procedure_2()
    function rectangular_pulse(t, T1)
        if t > -T1 && t < T1
            return 1.0
        else
            return 0.0
        end
    end

    fs = 1/100
    t = -10.0:fs:10.0
    w = -10.0:fs:10.0

    xt = rectangular_pulse.(t, 2)

    Xjw = fft(xt)
        
    function truncate_freq(Xjw, W, w)
        Xjw[w .> W] .= 0.0
        Xjw[w .< -W] .= 0.0
        return Xjw
    end

    Xjw_ = truncate_freq(Xjw, 5.0, w)

    x_t = ifft(Xjw_)

    et = xt - x_t

    #plot!(t, real(xt), label="x(t)")
    #plot!(t, real(x_t), label="x(t) truncated")

    function signal_energy(xt, W, Xjw, w)
        Xjw_ = truncate_freq(copy(Xjw), W, w)
        t = 0
        for i = Xjw
            if real(i) > 0.0
                t+=1
            end
        end
        println(W," ",t)
        x_t = ifft(Xjw_)

        et = xt - x_t

        return sum(et.^2)
    end

    W = 0.1:0.1:10
    eplot = signal_energy.(Ref(xt), W, Ref(Xjw), Ref(w))
    eplot
    #Plotting
    #plot!(W, real(eplot), label="Energy")
end

procedure_1()