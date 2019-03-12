classdef qpalm < handle
    % qpalm interface class for QPALM solver 
    % This class provides a complete interface to the C implementation
    % of the QPALM solver.
    %
    % qpalm Methods:
    %
    %   setup             - configure solver with problem data
    %   solve             - solve the QP
    %   warm_start        - set warm starting variables x and y TODO
    %
    %   default_settings  - create default settings structure
    %   current_settings  - get the current solver settings structure

    properties (SetAccess = private, Hidden = true)
%         objectHandle % Handle to underlying C instance
    end
   methods(Static) 
        %%
        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = qpalm_mex('default_settings');
        end
        
        %%
        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = qpalm_mex('constant', constant_name);
        end
        
    end
    methods
        %% Constructor - Create a new solver instance
        function this = qpalm(varargin)
            % Construct QPALM solver class
%             this.objectHandle = qpalm_mex('new', varargin{:});
        end

        %% Destructor - destroy the solver instance
        function delete(this)
            % Destroy QPALM solver class
            qpalm_mex('delete');
        end


        %%
        function varargout = setup(this, varargin)
            % SETUP configure solver with problem data
            %
            %   setup(Q,q,A,bmin,bmax,options)

            nargin = length(varargin);

            %dimension checks on user data. Mex function does not
            %perform any checks on inputs, so check everything here
            assert(nargin >= 5, 'incorrect number of inputs');
            [Q,q,A,bmin,bmax] = deal(varargin{1:5});

            %
            % Get problem dimensions
            %

            % Get number of variables n
            if (isempty(Q))
                if (~isempty(q))
                    n = length(q);
                else
                    if (~isempty(A))
                        n = size(A, 2);
                    else
                        error('The problem does not have any variables');
                    end
                end
            else
                n = size(Q, 1);
                assert(n==size(Q,2), 'Q must be a square matrix');
            end

            % Get number of constraints m
            if (isempty(A))
                m = 0;
            else
                m = size(A, 1);
            end

            %
            % Create sparse matrices and full vectors if they are empty
            %

            if (isempty(Q))
                Q = sparse(n, n);
%                 Q = zeros(n,n);
            else
                Q   = sparse(Q);
%                 Q = full(Q(:,:));
            end
            if (isempty(q))
                q = zeros(n, 1);
            else
                q   = full(q(:));
            end

            % Create proper constraints if they are not passed
            if (isempty(A) && (~isempty(bmin) || ~isempty(bmax))) || ...
                (~isempty(A) && (isempty(bmin) && isempty(bmax)))
                error('A must be supplied together with at least one bound l or u');
            end

            if (~isempty(A) && isempty(bmin))
                bmin = -Inf(m, 1);
            end

            if (~isempty(A) && isempty(bmax))
                bmax = Inf(m, 1);
            end

            if (isempty(A))
                A = sparse(m, n);
%                 A = zeros(m,n);
                bmin = -Inf(m, 1);
                bmax = Inf(m, 1);
            else
                bmin = full(bmin(:));
                bmax = full(bmax(:));
                A = sparse(A);
%                 A = full(A(:,:));
            end


            %
            % Check vector dimensions (not checked from the C solver)
            %

            assert(length(q)    == n, 'Incorrect dimension of q');
            assert(length(bmin) == m, 'Incorrect dimension of l');
            assert(length(bmax) == m, 'Incorrect dimension of u');

            %
            % Convert infinity values to QPALM_INFINITY
            %
            bmax = min(bmax, qpalm.constant('QPALM_INFTY'));
            bmin = max(bmin, -qpalm.constant('QPALM_INFTY'));


            %make a settings structure from the remainder of the arguments.
            %'true' means that this is a settings initialization, so all
            %parameter/values are allowed.  No extra inputs will result
            %in default settings being passed back
            theSettings = validateSettings(varargin{6:end});

            [varargout{1:nargout}] = qpalm_mex('setup',n,m,Q,q,A,bmin,bmax,theSettings);

        end


%         %%
% 
%         function warm_start(this, varargin)
%             % WARM_START warm start primal and/or dual variables
%             %
%             %   warm_start('x', x, 'y', y)
%             %
%             %   or warm_start('x', x)
%             %   or warm_start('y', y)
% 
% 
%             % Get problem dimensions
%             [n, m]  = get_dimensions(this);
% 
%             % Get data
%             allowedFields = {'x','y'};
% 
%             if(isempty(varargin))
%                 return;
%             elseif(length(varargin) == 1)
%                 if(~isstruct(varargin{1}))
%                     error('Single input should be a structure with new problem data');
%                 else
%                     newData = varargin{1};
%                 end
%             else % param / value style assumed
%                 newData = struct(varargin{:});
%             end
% 
%             %check for unknown fields
%             newFields = fieldnames(newData);
%             badFieldsIdx = find(~ismember(newFields,allowedFields));
%             if(~isempty(badFieldsIdx))
%                  error('Unrecognized input field ''%s'' detected',newFields{badFieldsIdx(1)});
%             end
% 
%             %get all of the terms.  Nonexistent fields will be passed
%             %as empty mxArrays
%             try x = double(full(newData.x(:))); catch x = []; end
%             try y = double(full(newData.y(:))); catch y = []; end
% 
%             % Check dimensions
%             assert(isempty(x) || length(x) == n, 'input ''x'' is the wrong size');
%             assert(isempty(y) || length(y) == m, 'input ''y'' is the wrong size');
% 
% 
%             % Decide which function to call
%             if (~isempty(x) && isempty(y))
%                 osqp_mex('warm_start_x', this.objectHandle, x);
%                 return;
%             end
% 
%             if (isempty(x) && ~isempty(y))
%                 osqp_mex('warm_start_y', this.objectHandle, y);
%             end
% 
%             if (~isempty(x) && ~isempty(y))
%                 osqp_mex('warm_start', this.objectHandle, x, y);
%             end
% 
%             if (isempty(x) && isempty(y))
%                 error('Unrecognized fields');
%             end
% 
%         end

        %%
        function varargout = solve(this, varargin)
            % SOLVE solve the QP

            nargoutchk(0,1);  %either return nothing (but still solve), or a single output structure
            [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = qpalm_mex('solve');
            if(nargout)
                varargout{1} = out;
            end
            delete(this);
            return;
        end

    end
end



function currentSettings = validateSettings(varargin)

%get the default settings
currentSettings = qpalm_mex('default_settings');

%no settings passed -> return defaults
if(isempty(varargin))
    return;
end

%check for structure style input
if(isstruct(varargin{1}))
    newSettings = varargin{1};
    assert(length(varargin) == 1, 'too many input arguments');
else
    newSettings = struct(varargin{:});
end

%get the qpalm settings fields
currentFields = fieldnames(currentSettings);

%get the requested fields in the update
newFields = fieldnames(newSettings);

%check for unknown parameters
badFieldsIdx = find(~ismember(newFields,currentFields));
if(~isempty(badFieldsIdx))
    error('Unrecognized solver setting ''%s'' detected',newFields{badFieldsIdx(1)});
end



%everything checks out - merge the newSettings into the current ones
for i = 1:length(newFields)
    currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
end


end


