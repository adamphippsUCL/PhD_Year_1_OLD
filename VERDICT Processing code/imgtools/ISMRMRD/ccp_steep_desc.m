function ccp_steep_desc
% Based on demo9
%
%
% GRAPPA reconstruction with the steepest descent step
% to illustrate the use of Acquisition Model projections.

if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

%try

    % define raw data source
    filenm = pref_uigetfile('ccp_steep_desc','filename') ;
   
    input_data = gadgetron.MR_Acquisitions(filenm);

    % pre-process acquisitions
    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoGadget'} ...
         {'RemoveROOversamplingGadget'}];
    acq_proc = gadgetron.AcquisitionsProcessor(prep_gadgets);
    fprintf('---\n pre-processing acquisitions...\n')
    preprocessed_data = acq_proc.process(input_data);
    pp_norm = preprocessed_data.norm();

    % perform reconstruction
    recon = gadgetron.MR_BasicGRAPPAReconstruction();
    recon.set_input(preprocessed_data);
    fprintf('---\n reconstructing...\n');
    recon.process();
    output = recon.get_output();
    complex_images = output.select(2);

    % compute coil sensitivity maps
    csms = gadgetron.MR_CoilSensitivityMaps();
    fprintf('---\n sorting acquisitions...\n')
    preprocessed_data.sort()
    fprintf('---\n calculating sensitivity maps...\n')
    csms.calculate(preprocessed_data)

    % create acquisition model based on the acquisition parameters
    % stored in preprocessed_data and image parameters stored in complex_images
    am = gadgetron.MR_AcquisitionModel(preprocessed_data, complex_images);
    am.set_coil_sensitivity_maps(csms)

    % use the acquisition model (forward projection) to simulate acquisitions
    fwd_data = am.forward(complex_images);
    fwd_norm = fwd_data.norm();
    % compute the difference between real and simulated acquisitions
    diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm);
    rr = diff.norm()/fwd_norm;
    fprintf('---\n reconstruction residual norm (rel): %e\n', rr)

    % try to improve the reconstruction by the steepest descent step
    g = am.backward(diff);
    w = am.forward(g);
    alpha = (g*g)/(w*w);
    r_complex_imgs = complex_images - g*alpha;

    % get real-valued reconstructed and refined images
    images = gadgetron.MR_extract_real_images(complex_images);
    r_imgs = gadgetron.MR_extract_real_images(r_complex_imgs);

    % plot images
    data1 = images.image_as_array(1);
    rdata1 = r_imgs.image_as_array(1);
    nim = images.number() ;
    data = zeros(size(data1,1), size(data1,2),nim) ;
    rdata = zeros(size(rdata1,1), size(rdata1,2),nim) ;
    for im = 1:nim
        data(:,:,im) = images.image_as_array(im);
        rdata(:,:,im) = r_imgs.image_as_array(im);
    end
        
    eshow(data,'Name','initial recon')
    eshow(rdata,'Name','refined recon')

    
    
% catch err
%     % display error information
%     fprintf('%s\n', err.message)
%     fprintf('error id is %s\n', err.identifier)
% end
