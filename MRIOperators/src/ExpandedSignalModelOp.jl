export ExpandedSignalModelOp

using SphericalHarmonicExpansions

function produ!(y::AbstractVector{T}, smaps::AbstractMatrix{T},higherOrderTrajectory::AbstractMatrix{Float32},sphericalHarmonicsBasis::AbstractMatrix{Float32}, x::AbstractVector{T}, numVox, numChan, numRead,numContr=1) where T

  x_ = reshape(x,numVox,numContr)
  y_ = reshape(y,numRead,numContr,numChan)

  @info("stop here produ ...")
  @infiltrate

  @time begin
    if !iszero(x_)
      for i_chan=1:numChan
        @inbounds @floop for i_read=1:numRead, i_vox = 1:numVox
            ϕ =  sum(higherOrderTrajectory[:,i_read[1]] .* sphericalHarmonicsBasis[i_vox[1],:])
            y_[i_read[1],numContr,i_chan] =  smaps[i_vox[1],i_chan] .* exp.(im.*(ϕ)).* x_[i_read[1],numContr] 
        end
      end
    end
  end

  @info("stop here produ ...")
  @infiltrate

  # x_ = x_ .* smaps 
  # for i=1:numChan
  #   y_[:,numContr,i] .= ϕ * x_[:,i]
  # end

  # @inbounds for i ∈ CartesianIndices(y_)
  #   y_[i] = x_[i[1],i[2]] * smaps[i[1],i[3]]
  # end

  return y
end

function ctprodu!(y::AbstractVector{T}, smapsC::AbstractMatrix{T},higherOrderTrajectory::AbstractMatrix{Float32},sphericalHarmonicsBasis::AbstractMatrix{Float32}, x::AbstractVector{T}, numVox, numChan, numRead, numContr=1) where T
  

  x_ = reshape(x,numRead,numContr,numChan)
  y_ = reshape(y,numVox,numContr)

  @info("stop here produ ...")
  @infiltrate

  @time begin
    if !iszero(x_)
      for i_chan=1:numChan
        @inbounds @floop for i_read=1:numRead, i_vox = 1:numVox
          ϕ =  conj(sum(higherOrderTrajectory[:,i_read[1]] .* sphericalHarmonicsBasis[i_vox[1],:]))
          y_[i_vox[1],numContr] =  smapsC[i_vox[1],i_chan] .* exp.(im.*(ϕ)).* x_[i_read[1],numContr,i_chan] 
        end
      end
    end
  end

  @info("stop here produ ...")
  @infiltrate

  # y_ .= 0

  # tmp = Array{ComplexF32,2}(undef,numVox,numChan)
  # for i=1:numChan
  #   # tmp[:,i] = transpose(smapsC[:,i]).*(transpose(x_[:,1,i]) * ϕ')
  #   tmp[:,i] = smapsC[:,i].*((adjoint(ϕ))*x_[:,1,i])
  #   # tmp[:,i] = smapsC[:,i].*(conj(ϕ)*x_[:,1,i])
  # end

  # y_ .= sum(tmp,dims=2)

  # @inbounds for i ∈ CartesianIndices(y_)
  #   y_[i] = x_[i[1],i[2]] * smaps[i[1],i[3]]
  # end


  return y
end


"""
ExpandedSignalModelOp(sensMaps::AbstractMatrix{T}, numContr=1)
(sensMaps::AbstractMatrix{T},OffResonanceMap::AbstractMatrix{T},shape::Tuple)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`
# Arguments
* `sensMaps`    - sensitivity maps ( 1. dim -> voxels, 2. dim-> coils)
* `numEchoes`   - number of contrasts to which the operator will be applied
"""
function ExpandedSignalModelOp(sensMaps::AbstractArray{T,4},OffResonanceMap::AbstractArray{T,3},higherOrderTrajectory::Matrix{Float32},shape::Tuple,FieldOfView::Tuple,times::Vector{Float32}, numContr=1) where T
    
    sensMaps = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
    OffResonanceMap = imag.(OffResonanceMap)
    OffResonanceMap = reshape(OffResonanceMap, div(length(OffResonanceMap),size(OffResonanceMap,4)),size(OffResonanceMap,4))

    numVox, numChan = size(sensMaps)
    numRead = size(times,1)
    sensMapsC = conj.(sensMaps) 

    # ToDo: Get variable to assign order l ..
    # ToDo: Not sure if I want to do the exponential here or in prod
    sphericalHarmonicsBasis = CalculateSphericalHarmonics(shape,FieldOfView;l=1)

    # # With OffResonanceMap
    # ToDo: Check units of OffResonanceMap
    # ϕ = exp.(im.*(OffResonanceMap.+sphericalHarmonicsBasis*higherOrderTrajectory))

    # Without OffResonanceMap
    # ϕ = exp.(im*Float32(3*pi).*(sphericalHarmonicsBasis*higherOrderTrajectory))
    # ϕ = exp.(im.*(transpose(Float32(2*pi).*higherOrderTrajectory)*transpose(sphericalHarmonicsBasis)))
    # ToDo: For some reson I need a factor of 10 more in sphericalHarmonicsBasis, do I need cm instead of mm?
    
    # @info("Stop...")
    # @infiltrate

    # ϕ = exp.(im.*(transpose(higherOrderTrajectory)*transpose(sphericalHarmonicsBasis.*Float32(1e1))))
    
    return LinearOperator{T}(numRead*numContr*numChan, numVox*numContr, false, false,
                         (res,x) -> produ!(res,sensMaps,higherOrderTrajectory,sphericalHarmonicsBasis,x,numVox,numChan,numRead,numContr),
                         nothing,
                         (res,x) -> ctprodu!(res,sensMapsC,higherOrderTrajectory,sphericalHarmonicsBasis,x,numVox,numChan,numRead,numContr))
end

# """
# ExpandedSignalModelOp(sensMaps::AbstractArray{T,4}, numContr=1)

# builds a `LinearOperator` which performs multiplication of a given image with
# the coil sensitivities specified in `sensMaps`
# # Arguments
# * `sensMaps`  - sensitivity maps ( 1.-3. dim -> voxels, 4. dim-> coils)
# * `numContr`  - number of contrasts to which the operator will be applied
# """
# function ExpandedSignalModelOp(sensMaps::AbstractArray{T,4},OffResonanceMap::AbstractArray{T,3},higherOrderTrajectory::Matrix{Float32},shape::Tuple,FieldOfView::Tuple,times::Vector{Float32}) where T #{T,D}
#   sensMaps = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
#   OffResonanceMap = reshape(OffResonanceMap, div(length(OffResonanceMap),size(OffResonanceMap,4)),size(OffResonanceMap,4))
#   return ExpandedSignalModelOp(sensMaps::AbstractMatrix{T},OffResonanceMap::AbstractMatrix{T},higherOrderTrajectory::Matrix{Float32},shape::Tuple,FieldOfView::Tuple,times::Vector{Float32})
# end

"""
CalculateSphericalHarmonics(params::Dict{Symbol,Any},SensitivityMap,OffResonanceMap,ks_traj; l::Int=3)

Creates Spherical harmonic basis functions
"""
function CalculateSphericalHarmonics(reconSize,FieldOfView; l::Int=3)

    lmax = (l+1).^2

    # SphericalHarmonicsBasis = zeros(reconSize[1],reconSize[2],reconSize[3],lmax)
    SphericalHarmonicsBasis = Array{Float32,4}(undef,reconSize[1],reconSize[2],reconSize[3],lmax)

    # Getting shperical harmonics basis functions
    @polyvar x y z

    l = 0:l
    m = -l[end]:l[end]

    # Creating cartesian coordinates
    (xx,yy,zz) = (LinRange.(-1,1,[reconSize[1],reconSize[2],reconSize[3]]).*FieldOfView./2)

    tmp = 1
    for  il=0:l[end]
        for im=-il:il
                f = rlylm(il,im,x,y,z)
                g = @fastfunc f
                for ix=eachindex(xx), iy=eachindex(yy), iz=eachindex(zz)
                    SphericalHarmonicsBasis[ix, iy, iz,tmp] = Base.invokelatest(g,xx[ix], yy[iy], zz[iz])
                end

                tmp = tmp+1
                # @info string("l=",il,", m=",im)
        end
    end

    SphericalHarmonicsBasis = reshape(SphericalHarmonicsBasis,:,lmax)

    return SphericalHarmonicsBasis
end

