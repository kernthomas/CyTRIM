#pragma once
#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>

template <typename T, std::size_t N>
class dirVector{
    private:
        std::array<T, N> v{};

    public:

        dirVector(std::initializer_list<T> list) {
            if (list.size() != N) {
                throw std::runtime_error("dirVector: wrong number of elements");
            }
            std::copy(list.begin(), list.end(), v.begin());
        }
        
        dirVector(){}

        T& operator[](std::size_t i){
            return v[i];
        }
        const T& operator[](std::size_t i) const{
            return v[i];
        }

        std::size_t size() const{
            return N;
        }

        dirVector operator+(const dirVector &w) const{

            dirVector<T, N> temp;

            if(w.size() != N){
                throw std::runtime_error("dirVector+: vector sizes don't match!");
            }
            
            for(int i=0; i<N; i++){
                temp[i] = v[i] + w[i];
            }
            return temp;
        }

        dirVector operator-(const dirVector &w) const{

            dirVector<T, N> temp;

            if(w.size() != N){
                throw std::runtime_error("dirVector+: vector sizes don't match!");
            }
            
            for(int i=0; i<N; i++){
                temp[i] = v[i] - w[i];
            }
            return temp;
        }

        dirVector operator*(const T factor) const{

            dirVector<T, N> result;

            for(int i=0; i<N; i++){
                result[i] = factor * v[i];
            }

            return result;
        }

        dirVector operator/(const T divisor) const{

            dirVector<T, N> result;

            for(int i=0; i<N; i++){
                result[i] = v[i] / divisor;
            }

            return result;
        }

        // Euclidean norm
        T norm() const {
            T result = 0;
            
            for(T element : v){
                result += (element * element);
            }

            result = sqrt(result);
            return result;
        }

        dirVector vectabs() const{
            dirVector<T, N> result;

            for(int i=0; i<N; i++){
                result[i] = std::abs(v[i]);
            }

            return result;
        }

        int argmin() const{
            auto it = std::min_element(v.begin(), v.end());

            if(it != v.end()){
                return static_cast<int>(std::distance(v.begin(), it));
            } else{
                return -1;
            }
        }
};