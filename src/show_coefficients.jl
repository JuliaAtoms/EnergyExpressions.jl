function showcoeff(io::IO, n::Number, show_sign::Bool, show_one::Bool=false)
    isimag(n) = !iszero(n) && isreal(im*n)
    isneg(n) = n < 0

    if isreal(n) && isneg(n) || isimag(n) && isneg(imag(n))
        write(io, "- ")
    else
        show_sign && write(io, "+ ")
    end
    if !(isone(n) || isone(-n))
        if !(isreal(n) || isimag(n))
            write(io, "($n)")
        elseif isimag(n)
            ii = imag(n)
            if !(isone(ii) || isone(-ii))
                write(io, "$(ii)im")
            else
                write(io, "im")
            end
        else
            write(io, string(n))
        end
    elseif show_one
        write(io, "1")
    end
end
