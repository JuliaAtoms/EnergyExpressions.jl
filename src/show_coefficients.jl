function showcoeff(io::IO, n::Number, show_sign::Bool, show_one::Bool=false)
    isimag(n) = !iszero(n) && isreal(im*n)
    isneg(n) = n < 0

    if isreal(n) && isneg(real(n)) || isimag(n) && isneg(imag(n))
        write(io, "- ")
    else
        show_sign && write(io, "+ ")
    end
    n = round(n, digits=6)
    if !(isone(n) || isone(-n))
        if !(isreal(n) || isimag(n))
            write(io, "($n)")
        elseif isimag(n)
            ii = imag(n)
            if !(isone(ii) || isone(-ii))
                write(io, "$(abs(ii))im")
            else
                write(io, "im")
            end
        else
            write(io, string(abs(n)))
        end
    elseif show_one
        write(io, "1")
    end
end

function showcoeff(io::IO, n, show_sign::Bool, _)
    show_sign && write(io, "+ ")
    print(io, n)
end
